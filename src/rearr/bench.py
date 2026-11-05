import os
import pathlib
import tempfile

import pandas as pd
import sh


def call_bash_func(bash_lib: str, *args):
    sh.bash(c=" ".join([".", bash_lib, "&&"] + args))


def bench_rearr(
    ref_file: os.PathLike,
    edit_file: os.PathLike,
    query_ratio: float,
    linux_file: os.PathLike,
    time_file: os.PathLike,
    alg_file: os.PathLike,
) -> None:
    ref_df = pd.read_csv(ref_file, header=0)
    edit_df = pd.read_csv(edit_file, header=0)
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)

        ref_df.assign(
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
        )[["rstart1", "ref1", "rend1", "rstart2", "ref2", "rend2"]].to_csv(
            td / "temp.ref", sep="\t", index=False, header=False
        )

        edit_df.assign(
            query=lambda df: df["mut1"] + df["random_insertion"] + df["mut2"],
            count=1,
        )[["query", "count", "ref_id"]].to_csv(
            td / "temp.query", sep="\t", index=False, header=False
        )

        call_bash_func(
            "src/rearr/bash_lib.sh",
            "bench_rearr",
            os.fspath(linux_file),
            os.fspath(time_file),
            os.fspath(td / "temp.ref"),
            os.fspath(td / "temp.quey"),
            os.fspath(alg_file),
        )
