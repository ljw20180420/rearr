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
from Bio.Seq import Seq


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


def parse_ref(ref_file: os.PathLike, query_ratio: float) -> tuple:
    ref_df = read_ref(ref_file, query_ratio)
    assert len(ref_df) == 1, "only single reference is supported"

    return (
        ref_df.loc[0, "ref1"],
        ref_df.loc[0, "ref2"],
        ref_df.loc[0, "rstart1"],
        ref_df.loc[0, "rend1"],
        ref_df.loc[0, "rstart2"],
        ref_df.loc[0, "rend2"],
    )


def parse_export(td: os.PathLike) -> tuple[float, float]:
    export_df = pd.read_csv(td / "export.csv", header=0)
    wall_time = export_df.loc[0, "time_wall_clock_mean"].item()
    memory = export_df.loc[0, "memory_peak_resident_mean"].item()

    return wall_time, memory


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

        wall_time, memory = parse_export(td)
        time_df = pd.read_csv(
            td / "temp.time", sep=": ", names=["name", "time"], engine="python"
        ).set_index("name")
        total_time = time_df.loc["total time", "time"].item()
        preprocess_time = time_df.loc["preprocess time", "time"].item()
        align_time = time_df.loc["align time", "time"].item()
        DP_time = time_df.loc["DP time", "time"].item()
        backtrack_time = time_df.loc["backtrack time", "time"].item()

    return (
        total_time,
        preprocess_time,
        align_time,
        DP_time,
        backtrack_time,
        wall_time,
        memory,
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
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        ref = ref1[:rend1] + ref2[rstart2:]
        sgRNA = ref[rend1 - 17 : rend1 + 3]
        cut1 = rend1
        cut2 = rstart2

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

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_amplicon(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        amplicon = (
            ref1[rstart1 - 6 : rend1].lower()
            + ref2[rstart2].upper()
            + ref2[rstart2 + 1 :].lower()
        )
        primer = amplicon[:6].upper()
        sgRNA = ref1[rend1 - 17 : rend1] + ref2[rstart2 : rstart2 + 3]
        cut1 = rend1 - rstart1 + 6
        cut2 = rstart2

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
                    src/rearr/amplicon.r \
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
                    rpos1=lambda df: df["rpos1"] - cut1 + rend1,
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

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_crisprvariants(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        ref = ref1[:rend1] + ref2[rstart2:]
        with open(td / "ref.fa", "w") as fd:
            fd.write(f">ref\n{ref}\n")
        cut1 = rend1
        cut2 = rstart2

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
            plus="+",
        )
        edit_df[["seqname", "query", "plus", "quality"]].stack().to_csv(
            td / "query.fq", index=False, header=False
        )

        sgRNA = ref[cut1 - 17 : cut1 + 3]
        bake_hyperfine(td, "crisprvariants")(
            textwrap.dedent(
                f"""\
                    bwa index -p {td / "ref"} {td / "ref.fa"}
			        bwa mem -T -9999 {td / "ref"} {td / "query.fq"} |
			        samtools view -bo {td / "query.bam"} -
			        src/rearr/crisprvariants.r \
                    {sgRNA} \
                    {td / "query.bam"} \
                    {cut1 - 17} \
                    > {td / "result"}"""
            ),
        )

        predict_df_slices = []
        with pysam.AlignmentFile(td / "query.bam") as bam:
            for read in bam:
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

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_ampliconDIVider(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        prefix = "AAAAAATGTAAAACGACGGCCAGT"
        suffix = "aagacac"
        ref = (prefix + ref1[rstart1:rend1] + ref2[rstart2:rend2] + suffix).upper()
        with open(td / "ref.fa", "w") as fd:
            fd.write(f">ref\n{ref}\n")
        cut1 = rend1 - rstart1 + len(prefix)
        cut2 = rstart2

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            padded_query=lambda df: df["query"].map(
                lambda query: f"{prefix}{query}{suffix}"
            ),
            rc_padded_query=lambda df: df["padded_query"].map(
                lambda padded_query: str(Seq(padded_query).reverse_complement()).upper()
            ),
            quality=lambda df: df["padded_query"].map(
                lambda padded_query: "~" * len(padded_query)
            ),
            plus="+",
        )
        edit_df[["seqname", "padded_query", "plus", "quality"]].stack().to_csv(
            td / "query.for.fq", index=False, header=False
        )
        edit_df[["seqname", "rc_padded_query", "plus", "quality"]].stack().to_csv(
            td / "query.rev.fq", index=False, header=False
        )

        with open(td / "query.bar", "w") as fd:
            fd.write(f"plate0\tA0\trefin_ref_{0}-{len(ref)}\tM\t0\tAAAAAA\n")

        ampliconDIVider_driver = pathlib.Path(
            "paper/tools/ampliconDIVider/ampliconDIVider_driver.sh"
        ).resolve()
        bake_hyperfine(td, "ampliconDIVider")(
            textwrap.dedent(
                f"""\
                    novoindex -k 6 -s 1 {td / "ref.nix"} {td / "ref.fa"}
			        novoalign -t 0,6 -g 10 -o SAM -d {td / "ref.nix"} -f {td / "query.for.fq"} {td / "query.rev.fq"} |
			        samtools view -o {td / "query.bam"}
			        cd {td}
			        {ampliconDIVider_driver} -b query.bar -f ref.fa -x ref.nix query.bam"""
            ),
        )

        predict_df_slices = []
        with pysam.AlignmentFile(td / "query.bam") as bam:
            for read in bam:
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
                        df["rpos1"] + reference_start - cut1 + rend1
                    ),
                )
                predict_df_slices.append(predict_df_slice)

        predict_df = pd.concat(predict_df_slices, ignore_index=True)[
            ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
        ].sort_values(by="idx", ascending=True, ignore_index=True)
        predict_df.to_csv(predict_file, index=False)

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_crispr_grant(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        ref = ref1[:rend1] + ref2[rstart2:]
        with open(td / "ref.fa", "w") as fd:
            fd.write(f">ref\n{ref}\n")
        cut1 = rend1
        cut2 = rstart2

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            rc_query=lambda df: df["query"].map(
                lambda query: str(Seq(query).reverse_complement())
            ),
            quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
            plus="+",
        )
        edit_df[["seqname", "query", "plus", "quality"]].stack().to_csv(
            td / "query.for.fq", index=False, header=False
        )
        edit_df[["seqname", "rc_query", "plus", "quality"]].stack().to_csv(
            td / "query.rev.fq", index=False, header=False
        )

        bake_hyperfine(td, "crispr_grant")(
            textwrap.dedent(
                f"""\
                    paper/tools/CRISPR-GRANT/bin/indel_analysis \
				        -1:{td / "query.for.fq"} \
				        -2:{td / "query.rev.fq"} \
				        -r:{td / "ref.fa"} \
				        -o:{td / "output"}"""
            ),
        )

        predict_df_slices = []
        with pysam.AlignmentFile(td / "output" / "4.Mapping_sorted.bam") as bam:
            for read in bam:
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

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_calculate_indel(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        sgRNA = ref1[rend1 - 17 : rend1] + ref2[rstart2 : rstart2 + 3]
        ref = ref1[rstart1:rend1] + ref2[rstart2:rend2]
        cut1 = rend1 - rstart1
        cut2 = rstart2

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
            plus="+",
        )
        edit_df[["seqname", "query", "plus", "quality"]].stack().to_csv(
            td / "query.fq", index=False, header=False
        )

        with open(td / "samplefile", "w") as fd:
            fd.write(f"query,{td / 'query.fq'},{sgRNA},{ref}\n")

        bake_hyperfine(td, "calculate_indel")(
            textwrap.dedent(
                f"""\
                    paper/tools/Screening_Protocols_manuscript/my_calculate_indel.py \
                        -i {td / "samplefile"} \
                        -no-m \
                        -o {td / "calc_indel_out.csv"} \
                        > {td / "result"}"""
            ),
        )

        predict_df_slices = []
        with open(td / "result", "r") as fd:
            next(fd)
            for line in fd:
                fields = line.rstrip().split()
                query = fields[0]
                ref_seq_coor = [int(s) for s in fields[1:]]
                predict_df_slice = coor2predict(
                    ref_seq_coor[: len(ref_seq_coor) // 2],
                    ref_seq_coor[len(ref_seq_coor) // 2 :],
                    cut1,
                    cut2,
                ).assign(
                    idx=int(query.removeprefix("seq")),
                    rpos1=lambda df: df["rpos1"] - cut1 + rend1,
                )
                predict_df_slices.append(predict_df_slice)

        predict_df = pd.concat(predict_df_slices, ignore_index=True)[
            ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
        ].sort_values(by="idx", ascending=True, ignore_index=True)
        predict_df.to_csv(predict_file, index=False)

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_SelfTarget(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        ref = ref1[:rend1] + ref2[rstart2:]
        with open(td / "oligo.fa", "w") as fd:
            fd.write(f">ref {rend1 + 3} FORWARD\n{ref}\n")
        cut1 = rend1
        cut2 = rstart2

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
            plus="+",
        )
        edit_df[["seqname", "query", "plus", "quality"]].stack().to_csv(
            td / "query.fq", index=False, header=False
        )

        bake_hyperfine(td, "SelfTarget")(
            textwrap.dedent(
                f"""\
                    paper/tools/SelfTarget-master/indel_analysis/indelmap/indelmap \
                        {td / "query.fq"} \
                        {td / "oligo.fa"} \
                        {td / "outputfile"} \
                        0"""
            ),
        )

        predict_df_slices = []
        with open(td / "outputfile", "r") as fd:
            next(fd)
            for line in fd:
                query, _, identifier, _ = line.split("\t", 3)
                if identifier == "-":
                    continue
                try:
                    (
                        _,
                        rpos1min,
                        rpos1max,
                        rpos2min,
                        rpos2max,
                        qpos1min,
                        qpos1max,
                        qpos2min,
                        qpos2max,
                    ) = identifier.split("|")
                except ValueError:
                    continue
                predict_df_slice = pd.DataFrame({
                    "idx": int(query.removeprefix("seq")),
                    "rpos1": [int(rpos1min), int(rpos1max)],
                    "rpos2": [int(rpos2min) - cut1 + cut2, int(rpos2max) - cut1 + cut2],
                    "qpos1": [int(qpos1min), int(qpos1max)],
                    "qpos2": [int(qpos2min), int(qpos2max)],
                })
                predict_df_slices.append(predict_df_slice)

        predict_df = pd.concat(predict_df_slices, ignore_index=True)[
            ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
        ].sort_values(by="idx", ascending=True, ignore_index=True)
        predict_df.to_csv(predict_file, index=False)

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_CRISPR_toolkit(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        user = "ljw"
        project = "rearr"
        barcode = ref1[rstart1 - 19 : rstart1]
        ref = barcode + ref1[rstart1:rend1] + ref2[rstart2:]
        (
            td / "Indel_searcher_2" / "Input" / user / "Reference" / project / "ref"
        ).mkdir(exist_ok=True, parents=True)
        with open(
            td
            / "Indel_searcher_2"
            / "Input"
            / user
            / "Reference"
            / project
            / "ref"
            / "Barcode.txt",
            "w",
        ) as fd:
            fd.write(f"{barcode}\n")
        with open(
            td
            / "Indel_searcher_2"
            / "Input"
            / user
            / "Reference"
            / project
            / "ref"
            / "Reference_sequence.txt",
            "w",
        ) as fd:
            fd.write(f"{ref}\n")
        with open(
            td
            / "Indel_searcher_2"
            / "Input"
            / user
            / "Reference"
            / project
            / "ref"
            / "Target_region.txt",
            "w",
        ) as fd:
            fd.write(f"{ref[len(barcode) :]}\n")
        cut1 = rend1 - rstart1 + len(barcode)
        cut2 = rstart2

        (td / "Indel_searcher_2" / "Input" / user / "FASTQ" / project / "query").mkdir(
            exist_ok=True, parents=True
        )
        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            barcode_query=lambda df: df["query"].map(lambda query: f"{barcode}{query}"),
            quality=lambda df: df["barcode_query"].map(
                lambda barcode_query: "~" * len(barcode_query)
            ),
            plus="+",
        )
        edit_df[["seqname", "barcode_query", "plus", "quality"]].stack().to_csv(
            td
            / "Indel_searcher_2"
            / "Input"
            / user
            / "FASTQ"
            / project
            / "query"
            / "query.fastq",
            index=False,
            header=False,
        )

        (td / "Indel_searcher_2" / "Core").mkdir(exist_ok=True, parents=True)
        sh.cp("paper/tools/CRISPR_toolkit/EDNAFULL", td)
        sh.cp(
            "paper/tools/CRISPR_toolkit/Indel_searcher_2/my_Run_indel_searcher.py",
            "paper/tools/CRISPR_toolkit/Indel_searcher_2/my_Indel_searcher_crispresso_hash.py",
            td / "Indel_searcher_2",
        )
        sh.cp(
            "paper/tools/CRISPR_toolkit/Core/__init__.py",
            "paper/tools/CRISPR_toolkit/Core/my_CoreSystem.py",
            td / "Indel_searcher_2" / "Core",
        )

        (td / "Indel_searcher_2" / "User" / user).mkdir(exist_ok=True, parents=True)
        with open(
            td / "Indel_searcher_2" / "User" / user / f"{project}.txt", "w"
        ) as fd:
            fd.write("query\tref\n")

        resso_python = (
            pathlib.Path(os.environ["HOME"])
            / ".pixi"
            / "envs"
            / "crispresso2"
            / "bin"
            / "python"
        )

        bake_hyperfine(td, "CRISPR_toolkit")(
            textwrap.dedent(
                f"""\
                    {resso_python} \
                    {td / "Indel_searcher_2" / "my_Run_indel_searcher.py"} \
				    --python {resso_python} \
				    --user {user} \
				    --project {project} \
				    --pam_type Cas9 \
				    --pam_pos Forward \
				    -t 15"""
            ),
        )

        predict_df_slices = []
        with open(td / "Indel_searcher_2" / "result", "r") as fd:
            for line in fd:
                seq_name, align_ref, align_seq = line.split("\t")
                refcoor, seqcoor = Alignment.parse_printed_alignment([
                    align_ref.encode(),
                    align_seq.encode(),
                ])[1]
                predict_df_slice = coor2predict(
                    refcoor, seqcoor, cut1 - len(barcode), cut2
                ).assign(
                    idx=int(seq_name.removeprefix("seq")),
                    rpos1=lambda df: df["rpos1"] + len(barcode) - cut1 + rend1,
                )
                predict_df_slices.append(predict_df_slice)

        predict_df = pd.concat(predict_df_slices, ignore_index=True)[
            ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
        ].sort_values(by="idx", ascending=True, ignore_index=True)
        predict_df.to_csv(predict_file, index=False)

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_CRISPR_Analytics(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        ref = ref1[:rend1] + ref2[rstart2:]
        sgRNA = ref[rend1 - 17 : rend1 + 3]
        with open(td / "ref.fa", "w") as fd:
            fd.write(f">ref\n{ref}\n")
        cut1 = rend1
        cut2 = rstart2

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
            plus="+",
        )
        edit_df[["seqname", "query", "plus", "quality"]].stack().to_csv(
            td / "query.fq", index=False, header=False
        )

        bake_hyperfine(td, "CRISPR_Analytics")(
            textwrap.dedent(
                f"""\
                    paper/tools/crispr-a_nextflow/bin/revComp-fasta.R \
                        {td / "ref.fa"} \
                        {td / "ref.correctOrient.fa"} \
                        {sgRNA}
                    minimap2 \
                        -d {td / "ref.correctOrient.mmi"} \
                        {td / "ref.correctOrient.fa"}
                    minimap2 -a -A 29 -B 17 -O 25 -E 2 \
                        {td / "ref.correctOrient.mmi"} \
                        {td / "query.fq"} \
                        > {td / "query.sam"}"""
            ),
        )

        predict_df_slices = []
        with pysam.AlignmentFile(td / "query.sam") as sam:
            for read in sam:
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

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_SIQ(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        ref = ref1[:rend1] + ref2[rstart2:]
        with open(td / "ref.fa", "w") as fd:
            fd.write(f">ref\n{ref}\n")
        left_flank = ref1[rstart1 - 15 : rstart1]
        right_flank = ref2[rend2 : rend2 + 15]
        cut1 = rend1
        cut2 = rstart2

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
            plus="+",
        )
        edit_df[["seqname", "query", "plus", "quality"]].stack().to_csv(
            td / "query.fq", index=False, header=False
        )

        bake_hyperfine(td, "SIQ")(
            textwrap.dedent(
                f"""\
                    java -classpath paper/tools/SIQ-1.8/SIQ/libraries/*:paper/tools/SIQ-1.8/SIQ/bin my.my_SIQ \
                        {td / "query.fq"} \
                        {td / "ref.fa"} \
                        {left_flank} \
                        {right_flank} \
                        {td / "output.tsv"}"""
            ),
        )

        predict_df_slices = []
        with open(td / "output.tsv", "r") as fd:
            for line in fd:
                seq_name, rpos1, rpos2, qpos1, qpos2 = line.split("\t")
                predict_df_slice = pd.DataFrame({
                    "idx": [int(seq_name.removeprefix("seq"))],
                    "rpos1": [int(rpos1)],
                    "rpos2": [int(rpos2) - cut1 + cut2],
                    "qpos1": [int(qpos1)],
                    "qpos2": [int(qpos2)],
                })
                predict_df_slices.append(predict_df_slice)

        predict_df = pd.concat(predict_df_slices, ignore_index=True)[
            ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
        ].sort_values(by="idx", ascending=True, ignore_index=True)
        predict_df.to_csv(predict_file, index=False)

        wall_time, memory = parse_export(td)

    return wall_time, memory


def bench_PEM_Q(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float, float]:
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(ref_file, query_ratio)
        ref = ref1[:rend1] + ref2[rstart2:]
        with open(td / "ref.fa", "w") as fd:
            fd.write(f">ref\n{ref}\n")
        cut1 = rend1
        cut2 = rstart2

        edit_df = read_edit(edit_file).assign(
            seqname=lambda df: df.index.map(lambda index: f"@seq{index}"),
            quality=lambda df: df["query"].map(lambda query: "~" * len(query)),
            plus="+",
        )
        edit_df[["seqname", "query", "plus", "quality"]].stack().to_csv(
            td / "query.fq", index=False, header=False
        )

        (td / "ref").mkdir(exist_ok=True, parents=True)
        sh.bwa.index("-p", td / "ref" / "ref", td / "ref.fa")
        bake_hyperfine(td, "PEM_Q")(
            textwrap.dedent(
                f"""\
                    BWA_INDEX={td} SAM_DIR={td} python paper/tools/PEM-Q-1.0.0/main/align_make_v5.1.py ref {td / "query.fq"}
                    SAM_DIR={td} python paper/tools/PEM-Q-1.0.0/main/define_indel_v5.1_mpf.py "_s" {cut1} +"""
            ),
        )

        idxs = []
        rpos1s = []
        rpos2s = []
        qpos1s = []
        qpos2s = []

        index_df = pd.read_csv(
            td / "raw" / "_s_indel_cutoff.tab",
            sep="\t",
            header=0,
            keep_default_na=False,
        )
        for (
            qname,
            rpos1,
            rpos2,
            strand,
            qpos1,
            qpos2,
            microhomology,
        ) in index_df[
            [
                "Qname",
                "Bait_end",
                "Prey_start",
                "Strand",
                "B_Qend",
                "Qstart",
                "Microhomolog",
            ]
        ].itertuples(index=False):
            if strand == "-":
                continue
            idx = int(qname.removeprefix("seq"))
            rpos2 = rpos2 - cut1 + cut2
            rpos1 -= len(microhomology)
            qpos1 -= len(microhomology)
            for shift in range(len(microhomology) + 1):
                idxs.append(idx)
                rpos1s.append(rpos1 + shift)
                rpos2s.append(rpos2 + shift)
                qpos1s.append(qpos1 + shift)
                qpos2s.append(qpos2 + shift)

        germline_addup_df = pd.read_csv(
            td / "raw" / "_s_Germline_addup_indel.tab",
            sep="\t",
            header=0,
            keep_default_na=False,
        )
        for (
            qname,
            rpos1,
            rpos2,
            strand,
            qpos1,
            qpos2,
            microhomology,
        ) in germline_addup_df[
            [
                "Qname",
                "Bait_end",
                "Prey_start",
                "Strand",
                "B_Qend",
                "Qstart",
                "Microhomolog",
            ]
        ].itertuples(index=False):
            if strand == "-":
                continue
            idx = int(qname.removeprefix("seq"))
            rpos2 = rpos2 - cut1 + cut2
            rpos1 -= len(microhomology)
            qpos1 -= len(microhomology)
            qpos1 += cut1 - rpos1
            qpos2 += cut2 - rpos2
            rpos1 = cut1
            rpos2 = cut2
            idxs.append(idx)
            rpos1s.append(rpos1 + shift)
            rpos2s.append(rpos2 + shift)
            qpos1s.append(qpos1 + shift)
            qpos2s.append(qpos2 + shift)

        predict_df_slices = [
            pd.DataFrame({
                "idx": idxs,
                "rpos1": rpos1s,
                "rpos2": rpos2s,
                "qpos1": qpos1s,
                "qpos2": qpos2s,
            })
        ]
        for pem_q_out_file in ["_s_Germline.tab", "_s_Substitution.tab"]:
            pem_q_df = pd.read_csv(td / "raw" / pem_q_out_file, sep="\t", header=0)
            for qname, cigar, reference_start in zip(
                pem_q_df["Qname"], pem_q_df["Cigar"], pem_q_df["Position"]
            ):
                reference_start -= 1
                read = pysam.AlignedSegment()
                read.cigarstring = cigar
                refcoor, seqcoor = sam2coor(read)
                predict_df_slice = coor2predict(
                    refcoor, seqcoor, cut1 - reference_start, cut2
                ).assign(
                    idx=int(qname.removeprefix("seq")),
                    rpos1=lambda df, reference_start=reference_start: (
                        df["rpos1"] + reference_start
                    ),
                )

                predict_df_slices.append(predict_df_slice)

        predict_df = pd.concat(predict_df_slices, ignore_index=True)[
            ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
        ].sort_values(by="idx", ascending=True, ignore_index=True)
        predict_df.to_csv(predict_file, index=False)

        wall_time, memory = parse_export(td)

    return wall_time, memory


def get_bench_func(lib: str) -> Callable:
    bench_funcs = {
        "rearr": bench_rearr,
        "resso": bench_resso,
        "amplicon": bench_amplicon,
        "crisprvariants": bench_crisprvariants,
        "ampliconDIVider": bench_ampliconDIVider,
        "crispr_grant": bench_crispr_grant,
        "calculate_indel": bench_calculate_indel,
        "SelfTarget": bench_SelfTarget,
        "CRISPR_toolkit": bench_CRISPR_toolkit,
        "CRISPR-Analytics": bench_CRISPR_Analytics,
        "SIQ": bench_SIQ,
        "PEM-Q": bench_PEM_Q,
    }

    return bench_funcs[lib]
