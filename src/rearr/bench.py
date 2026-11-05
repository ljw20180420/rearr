import os
import pathlib
import re
import tempfile

import Bio
import pandas as pd
import pysam
import sh


def call_bash_func(bash_lib: str, *args):
    sh.bash(c=" ".join([".", bash_lib, "&&"] + [*args]))


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


def cigar2coor(cigar):
    tuplestring = "MIDNSHP=XB"
    refcoor, seqcoor = [0], [0]
    pAS = pysam.AlignedSegment()
    pAS.cigarstring = cigar
    for op, sz in pAS.cigartuples:
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


def bench_rearr(
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
        call_bash_func(
            "src/rearr/bash_lib.sh",
            "bench_rearr",
            os.fspath(td / "temp.linux"),
            os.fspath(td / "temp.time"),
            os.fspath(td / "temp.ref"),
            os.fspath(td / "temp.query"),
            os.fspath(td / "temp.alg"),
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

        linux_df = pd.read_csv(
            td / "temp.linux", names=["real time", "user time", "system time", "RSS"]
        )
        time_df = pd.read_csv(
            td / "temp.time", sep=": ", names=["name", "time"], engine="python"
        ).set_index("name")

        return (
            time_df.loc["total time", "time"].item(),
            time_df.loc["preprocess time", "time"].item(),
            time_df.loc["align time", "time"].item(),
            time_df.loc["DP time", "time"].item(),
            time_df.loc["backtrack time", "time"].item(),
            linux_df.loc[0, "real time"].item(),
            linux_df.loc[0, "user time"].item(),
            linux_df.loc[0, "system time"].item(),
            linux_df.loc[0, "RSS"].item(),
        )


def bench_resso(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    predict_file: os.PathLike,
) -> tuple[float]:
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
        call_bash_func(
            "src/rearr/bash_lib.sh",
            "bench_resso",
            os.fspath(td / "temp.linux"),
            ref,
            sgRNA,
            os.fspath(td / "temp.query"),
            os.fspath(td / "output"),
        )
        zip_file = sh.find(td / "output", "-name", "*.zip")
        sh.unzip("-q", "-o", zip_file, d=td / "output")
        predict_df_slices = []
        with open(td / "output" / "Alleles_frequency_table.txt", "r") as fd:
            for line in fd:
                align_seq, align_ref, _ = line.split("\t", 2)
                refcoor, seqcoor = Bio.Align.Alignment.parse_printed_alignment([
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
                right=edit_df["seqname", "query"],
                on="query",
                how="left",
                validate="many_to_one",
            )
            .assign(idx=lambda df: df["seqname"].str.removeprefix("@seq").astype(int))[
                ["idx", "rpos1", "rpos2", "qpos1", "qpos2"]
            ]
            .sort_values(by="idx", ascending=True)
        )
        predict_df.to_csv(predict_file, index=False)

        linux_df = pd.read_csv(
            td / "temp.linux", names=["real time", "user time", "system time", "RSS"]
        )

        return (
            linux_df.loc[0, "real time"].item(),
            linux_df.loc[0, "user time"].item(),
            linux_df.loc[0, "system time"].item(),
            linux_df.loc[0, "RSS"].item(),
        )
