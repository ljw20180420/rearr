import itertools
import os
import pathlib
import re
import tempfile
import textwrap
from collections.abc import Callable

import pandas as pd
import pysam
import sh
from Bio.Align import Alignment


def read_ref(ref_file: os.PathLike, query_ratio: float) -> pd.DataFrame:
    ref_df = pd.read_csv(ref_file, header=0)
    ref_df = ref_df.assign(
        rstart1=lambda df: df["ref1"].map(
            lambda ref1: int(len(ref1) * (0.5 - query_ratio / 2))
        ),
        rend1=lambda df: df["ref1"].map(
            lambda ref1: int(len(ref1) * (0.5 + query_ratio / 2))
        ),
        rstart2=lambda df: df["ref2"].map(
            lambda ref2: int(len(ref2) * (0.5 - query_ratio / 2))
        ),
        rend2=lambda df: df["ref2"].map(
            lambda ref2: int(len(ref2) * (0.5 + query_ratio / 2))
        ),
    )[["rstart1", "ref1", "rend1", "rstart2", "ref2", "rend2"]]

    return ref_df


def read_edit(edit_file: os.PathLike) -> pd.DataFrame:
    edit_df = pd.read_csv(edit_file, header=0, keep_default_na=False).assign(
        query=lambda df: df["mut1"] + df["random_insertion"] + df["mut2"],
    )

    return edit_df


def coor2predict(refcoor, seqcoor, cut1, cut2) -> pd.DataFrame:
    segs = []
    for i in range(len(refcoor) - 1):
        if refcoor[i + 1] > refcoor[i] and seqcoor[i + 1] > seqcoor[i]:
            segs.append([refcoor[i], refcoor[i + 1], seqcoor[i], seqcoor[i + 1]])

    rpos1s = []
    rpos2s = []
    qpos1s = []
    qpos2s = []
    for i in range(len(segs) - 1):
        rpos1s.append(segs[i][1])
        rpos2s.append(segs[i + 1][0] - cut1 + cut2)
        qpos1s.append(segs[i][3])
        qpos2s.append(segs[i + 1][2])
    if not rpos1s:
        seqpos = segs[0][2] + cut1 - segs[0][0]
        rpos1s.append(cut1)
        rpos2s.append(cut2)
        qpos1s.append(seqpos)
        qpos2s.append(seqpos)

    return pd.DataFrame({
        "rpos1": rpos1s,
        "rpos2": rpos2s,
        "qpos1": qpos1s,
        "qpos2": qpos2s,
    })


def sam2coor(read: pysam.AlignedSegment) -> tuple[list[int], list[int]]:
    tuplestring = "MIDNSHP=XB"
    refcoor, seqcoor = [0], [0]
    for op, sz in read.cigartuples:
        refnc = refcoor[-1] + sz if tuplestring[op] in "M=XDN" else refcoor[-1]
        seqnc = seqcoor[-1] + sz if tuplestring[op] in "M=XISH" else seqcoor[-1]
        refcoor.append(refnc)
        seqcoor.append(seqnc)

    return refcoor, seqcoor


class RearrRQposParser:
    def __init__(self) -> None:
        self.lower_pattern = re.compile(r"[acgtn]")
        self.upper_pattern = re.compile(r"[ACGTN]")
        self.mid1_pattern = re.compile(r"[acgtn]-*[acgtn]")
        self.qend1_pattern = re.compile(r"[ACGTN]-*$")

    def __call__(self, refline: str, queryline: str) -> tuple[int, int, int, int]:
        upos = self.lower_pattern.search(refline).start()
        mid1 = self.mid1_pattern.search(refline, pos=upos + 1).start() + 1
        qstart1 = self.upper_pattern.search(queryline, pos=upos).start()
        qend1 = self.qend1_pattern.search(queryline[:mid1], pos=qstart1 + 1).start() + 1
        mid2 = self.lower_pattern.search(refline, pos=mid1).start()
        qstart2 = self.upper_pattern.search(queryline, pos=mid2).start()
        rpos1 = len(refline[:qend1].replace("-", ""))
        qpos1 = len(queryline[:qend1].replace("-", ""))
        rpos2 = len(refline[mid2:qstart2].replace("-", ""))
        qpos2 = len(queryline[:qstart2].replace("-", ""))

        return rpos1, rpos2, qpos1, qpos2


def bake_hyperfine(td: os.PathLike, command_name: str):
    td = pathlib.Path(os.fspath(td))

    return sh.hyperfine.bake(
        "--metrics",
        "time_wall_clock:s,memory_peak_resident:kB",
        "--runs",
        1,
        "--shell",
        "/bin/bash",
        "--export-csv",
        td / "export.csv",
        "--command-name",
        command_name,
    )


def bench_rearr_full(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        read_ref(ref_file, query_ratio).to_csv(
            td / "temp.ref", sep="\t", index=False, header=False
        )
        read_edit(edit_file).assign(
            count=1,
        )[["query", "count", "ref_id"]].to_csv(
            td / "temp.query", sep="\t", index=False, header=False
        )
        bake_hyperfine(td, "rearr")(
            textwrap.dedent(
                f"""\
                    core/Rearrangement/build/rearrangement \
                        <{td / "temp.query"} \
                        3<{td / "temp.ref"} \
                        >{td / "temp.alg"} \
                        4>{td / "temp.time"}"""
            ),
        )

        alg_df = pd.DataFrame(
            pd.read_csv(td / "temp.alg", header=None).values.reshape(-1, 3),
            columns=["header", "refline", "queryline"],
        )
        alg_df = alg_df.assign(
            idx=lambda df: (
                df["header"].str.split("\t", expand=True)[0].astype({0: int}) - 1
            )
        )
        rearr_rqpos_parser = RearrRQposParser()
        rpos1s = []
        rpos2s = []
        qpos1s = []
        qpos2s = []
        for refline, queryline in zip(alg_df["refline"], alg_df["queryline"]):
            rpos1, rpos2, qpos1, qpos2 = rearr_rqpos_parser(refline, queryline)
            rpos1s.append(rpos1)
            rpos2s.append(rpos2)
            qpos1s.append(qpos1)
            qpos2s.append(qpos2)

        alg_df.assign(
            rpos1=rpos1s,
            rpos2=rpos2s,
            qpos1=qpos1s,
            qpos2=qpos2s,
        )[["idx", "rpos1", "rpos2", "qpos1", "qpos2"]].to_csv(predict_file, index=False)

        export_df = pd.read_csv(td / "export.csv", header=0)
        time_df = pd.read_csv(
            td / "temp.time", sep=": ", names=["name", "time"], engine="python"
        ).set_index("name")

        return (
            time_df.loc["total time", "time"].item(),
            time_df.loc["preprocess time", "time"].item(),
            time_df.loc["align time", "time"].item(),
            time_df.loc["DP time", "time"].item(),
            time_df.loc["backtrack time", "time"].item(),
            export_df.loc[0, "time_wall_clock_mean"].item(),
            export_df.loc[0, "memory_peak_resident_mean"].item(),
        )


def bench_rearr(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    (
        _,
        _,
        _,
        _,
        _,
        wall_time,
        memory,
    ) = bench_rearr_full(ref_file, edit_file, query_ratio, predict_file)

    return wall_time, memory


def bench_resso(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref_df = read_ref(ref_file, query_ratio)
        assert len(ref_df) == 1, "only single reference is supported"
        ref1 = ref_df.loc[0, "ref1"]
        cut1 = ref_df.loc[0, "rend1"]
        ref2 = ref_df.loc[0, "ref2"]
        cut2 = ref_df.loc[0, "rstart2"]
        ref = ref1[:cut1] + ref2[cut2:]
        sgRNA = ref[cut1 - 17 : cut1 + 3]

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
            plus="+",
        )
        edit_df[["seqname", "query", "plus", "quality"]].stack().to_csv(
            td / "temp.query", index=False, header=False
        )
        (td / "output").mkdir(exist_ok=True, parents=True)
        bake_hyperfine(td, "resso")(
            textwrap.dedent(
                f"""\
                    CRISPResso \
		            -r1 {td / "temp.query"} \
		            -a {ref} \
		            -amas 0 \
		            --default_min_aln_score 0 \
		            -g {sgRNA} \
		            -o {td / "output"} \
		            --suppress_report \
		            --suppress_plots"""
            ),
        )
        zip_file = sh.find(td / "output", "-name", "*.zip").rstrip("\n")
        sh.unzip("-q", "-o", zip_file, d=td / "output")
        predict_df_slices = []
        with open(td / "output" / "Alleles_frequency_table.txt", "r") as fd:
            next(fd)
            for line in fd:
                align_seq, align_ref, _ = line.split("\t", 2)
                refcoor, seqcoor = Alignment.parse_printed_alignment([
                    align_ref.encode(),
                    align_seq.encode(),
                ])[1]
                predict_df_slice = coor2predict(refcoor, seqcoor, cut1, cut2)
                predict_df_slice = predict_df_slice.assign(
                    query=align_seq.replace("-", "")
                )
                predict_df_slices.append(predict_df_slice)

        predict_df = (
            pd
            .concat(predict_df_slices, ignore_index=True)
            .merge(
                right=edit_df[["seqname", "query"]],
                on="query",
                how="left",
                validate="many_to_one",
            )
            .assign(idx=lambda df: df["seqname"].str.removeprefix("@seq").astype(int))[
                ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
            ]
            .sort_values(by="idx", ascending=True, ignore_index=True)
        )
        predict_df.to_csv(predict_file, index=False)

        export_df = pd.read_csv(td / "export.csv", header=0)
        wall_time = export_df.loc[0, "time_wall_clock_mean"].item()
        memory = export_df.loc[0, "memory_peak_resident_mean"].item()

        return wall_time, memory


def bench_amplicon(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref_df = read_ref(ref_file, query_ratio)
        assert len(ref_df) == 1, "only single reference is supported"
        ref1 = ref_df.loc[0, "ref1"]
        ref2 = ref_df.loc[0, "ref2"]
        rstart1 = ref_df.loc[0, "rstart1"]
        rend1 = ref_df.loc[0, "rend1"]
        rstart2 = ref_df.loc[0, "rstart2"]

        amplicon = (
            ref1[rstart1 - 6 : rend1].lower()
            + ref2[rstart2].upper()
            + ref2[rstart2 + 1 :].lower()
        )
        primer = amplicon[:6].upper()
        sgRNA = ref1[rend1 - 17 : rend1] + ref2[rstart2 : rstart2 + 3]

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            primer_query=lambda df: df["query"].map(lambda query: f"{primer}{query}"),
            quality=lambda df: df["primer_query"].map(
                lambda prime_query: "~" * len(prime_query)
            ),
            plus="+",
        )
        (td / "fastq").mkdir(exist_ok=True, parents=True)
        edit_df[["seqname", "primer_query", "plus", "quality"]].stack().to_csv(
            td / "fastq" / "query.fq", index=False, header=False
        )

        with open(td / "config.csv", "w") as fd:
            fd.write(
                f"ID,Barcode,Forward_Reads,Reverse_Reads,Group,Control,guideRNA,Forward_Primer,Reverse_Primer,Direction,Amplicon,Donor\nID_1,barcode_1,query.fq,,group_1,0,{sgRNA},{primer},,0,{amplicon},\n",
            )

        bake_hyperfine(td, "amplicon")(
            textwrap.dedent(
                f"""\
                    src/rearr/AMP.r \
				    {td / "config.csv"} \
				    {td / "fastq"} \
                    {td / "alignment.txt"}"""
            ),
        )

        unique_seq_query_df = (
            edit_df[["seqname", "primer_query"]]
            .groupby("primer_query", as_index=False)
            .agg(
                count=pd.NamedAgg(column="primer_query", aggfunc="size"),
                seqname=pd.NamedAgg(column="seqname", aggfunc="first"),
            )
            .sort_values(
                by=["count", "primer_query"], ascending=[False, True], ignore_index=True
            )
        )
        amplicon_start = rstart1 - 6
        cut1 = rend1 - amplicon_start
        cut2 = rstart2
        predict_df_slices = []
        with open(td / "alignment.txt", "r") as fd:
            for (_, align_ref, align_seq, _), primer_query in zip(
                itertools.batched(fd, 4), unique_seq_query_df["primer_query"]
            ):
                refcoor, seqcoor = Alignment.parse_printed_alignment([
                    align_ref.encode(),
                    align_seq.encode(),
                ])[1]
                predict_df_slice = coor2predict(refcoor, seqcoor, cut1, cut2)
                predict_df_slice = predict_df_slice.assign(
                    rpos1=lambda df: df["rpos1"] + amplicon_start,
                    qpos1=lambda df: df["qpos1"] - len(primer),
                    qpos2=lambda df: df["qpos2"] - len(primer),
                    primer_query=primer_query,
                )
                predict_df_slices.append(predict_df_slice)

        predict_df = (
            pd
            .concat(predict_df_slices, ignore_index=True)
            .merge(
                right=unique_seq_query_df[["seqname", "primer_query"]],
                on="primer_query",
                how="left",
                validate="many_to_one",
            )
            .assign(idx=lambda df: df["seqname"].str.removeprefix("@seq").astype(int))[
                ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
            ]
            .sort_values(by="idx", ascending=True, ignore_index=True)
        )
        predict_df.to_csv(predict_file, index=False)

        export_df = pd.read_csv(td / "export.csv", header=0)
        wall_time = export_df.loc[0, "time_wall_clock_mean"].item()
        memory = export_df.loc[0, "memory_peak_resident_mean"].item()

        return wall_time, memory


def bench_crisprvariants(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref_df = read_ref(ref_file, query_ratio)
        assert len(ref_df) == 1, "only single reference is supported"
        ref1 = ref_df.loc[0, "ref1"]
        ref2 = ref_df.loc[0, "ref2"]
        cut1 = ref_df.loc[0, "rend1"]
        cut2 = ref_df.loc[0, "rstart2"]
        ref = ref1[:cut1] + ref2[cut2:]
        with open(td / "ref.fa", "w") as fd:
            fd.write(f">ref\n{ref}\n")

        edit_df = (
            read_edit(edit_file)
            .assign(
                seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
                quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
                plus="+",
            )[["seqname", "query", "plus", "quality"]]
            .stack()
            .to_csv(td / "query.fq", index=False, header=False)
        )

        sgRNA = ref[cut1 - 17 : cut1 + 3]
        bake_hyperfine(td, "amplicon")(
            textwrap.dedent(
                f"""\
                    bwa index -p {td / "ref"} {td / "ref.fa"}
			        bwa mem -T -9999 {td / "ref"} {td / "query.fq"} |
			        samtools sort -o {td / "query.bam"}
			        samtools index -b {td / "query.bam"}
			        src/rearr/crisprvariants.r \
                    {sgRNA} \
                    {td / "query.bam"} \
                    {cut1 - 17} \
                    > {td / "result"}"""
            ),
        )

        predict_df_slices = []
        with pysam.AlignmentFile(td / "query.bam") as bam:
            for read in bam.fetch():
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_reverse
                ):
                    continue

                refcoor, seqcoor = sam2coor(read)
                predict_df_slice = coor2predict(
                    refcoor, seqcoor, cut1 - read.reference_start, cut2
                ).assign(
                    idx=int(read.query_name.removeprefix("seq")),
                    rpos1=lambda df, reference_start=read.reference_start: (
                        df["rpos1"] + reference_start
                    ),
                )
                predict_df_slices.append(predict_df_slice)

        predict_df = pd.concat(predict_df_slices, ignore_index=True)[
            ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
        ].sort_values(by="idx", ascending=True, ignore_index=True)
        predict_df.to_csv(predict_file, index=False)

        export_df = pd.read_csv(td / "export.csv", header=0)
        wall_time = export_df.loc[0, "time_wall_clock_mean"].item()
        memory = export_df.loc[0, "memory_peak_resident_mean"].item()

        return wall_time, memory


def get_bench_func(lib: str) -> Callable:
    bench_funcs = {
        "rearr": bench_rearr,
        "resso": bench_resso,
        "amplicon": bench_amplicon,
        "crisprvariants": bench_crisprvariants,
    }

    return bench_funcs[lib]
