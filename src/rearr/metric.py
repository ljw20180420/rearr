import os
import pathlib

import pandas as pd

from .bench import parse_ref, read_edit, read_ref


def expand_to_mh_equivalence(
    predict_df: pd.DataFrame, ref1: str, ref2: str, cut1: int, cut2: int, lib: str
) -> pd.DataFrame:
    rpos2_correct = 0
    if lib != "rearr":
        ref = ref1[:cut1] + ref2[cut2:]
        ref1 = ref
        ref2 = ref
        rpos2_correct = -cut2 + cut1

    idxs = []
    rpos1s = []
    rpos2s = []
    qpos1s = []
    qpos2s = []
    for idx, rpos1, rpos2, qpos1, qpos2 in predict_df.itertuples(index=False):
        lext = 0
        rext = 0
        if qpos1 == qpos2:
            while True:
                left_rpos1 = rpos1 - lext - 1
                left_rpos2 = rpos2 + rpos2_correct - lext - 1
                if (
                    left_rpos1 < 0
                    or left_rpos2 < 0
                    or ref1[left_rpos1] != ref2[left_rpos2]
                ):
                    break
                lext += 1
            while True:
                right_rpos1 = rpos1 + rext
                right_rpos2 = rpos2 + rpos2_correct + rext
                if (
                    right_rpos1 >= len(ref1)
                    or right_rpos2 >= len(ref2)
                    or ref1[right_rpos1] != ref2[right_rpos2]
                ):
                    break
                rext += 1

        for shift in range(-lext, rext + 1):
            idxs.append(idx)
            rpos1s.append(rpos1 + shift)
            rpos2s.append(rpos2 + shift)
            qpos1s.append(qpos1 + shift)
            qpos2s.append(qpos2 + shift)

    return pd.DataFrame({
        "idx": idxs,
        "rpos1": rpos1s,
        "rpos2": rpos2s,
        "qpos1": qpos1s,
        "qpos2": qpos2s,
    })


def select_best_indel(
    predict_file: os.PathLike,
    ref_file: os.PathLike,
    query_ratio: float,
    edit_file: os.PathLike,
    select_file: os.PathLike,
    lib: str,
) -> None:
    ref1, ref2, rstart1, rend1, rstart2, rend2 = parse_ref(
        read_ref(ref_file, query_ratio)
    )
    predict_df = pd.read_csv(predict_file, header=0)
    expand_df = expand_to_mh_equivalence(
        predict_df=predict_df,
        ref1=ref1,
        ref2=ref2,
        cut1=rend1,
        cut2=rstart2,
        lib=lib,
    )
    edit_df = pd.read_csv(edit_file, header=0, keep_default_na=False)

    edit_df = edit_df.assign(
        qpos1=lambda df: df["mut1"].str.len(),
        qpos2=lambda df: df["qpos1"] + df["random_insertion"].str.len(),
    ).reset_index(names=["idx"])

    expand_df = expand_df.merge(
        right=edit_df,
        how="left",
        on="idx",
        suffixes=("_p", "_e"),
        validate="many_to_one",
    )
    expand_df = expand_df.assign(
        total_err=lambda df: (
            (df["rpos1_p"] - df["rpos1_e"]).abs()
            + (df["rpos2_p"] - df["rpos2_e"]).abs()
            + (df["qpos1_p"] - df["qpos1_e"]).abs()
            + (df["qpos2_p"] - df["qpos2_e"]).abs()
        ),
        min_total_err=lambda df: df.groupby("idx")["total_err"].transform("min"),
    )
    expand_df = expand_df.query("total_err == min_total_err").reset_index(drop=True)
    expand_df = expand_df.assign(
        center_err=lambda df: (
            df["rpos1_p"]
            + df["rpos2_p"]
            + df["qpos1_p"]
            + df["qpos2_p"]
            - df["rpos1_e"]
            - df["rpos2_e"]
            - df["qpos1_e"]
            - df["qpos2_e"]
        ).abs(),
        min_center_err=lambda df: df.groupby("idx")["center_err"].transform("min"),
    )
    expand_df = expand_df.query("center_err == min_center_err").reset_index(drop=True)
    expand_df = expand_df.groupby("idx", as_index=False).first()
    expand_df[
        [
            "idx",
            "rpos1_p",
            "rpos2_p",
            "qpos1_p",
            "qpos2_p",
            "rpos1_e",
            "rpos2_e",
            "qpos1_e",
            "qpos2_e",
        ]
    ].to_csv(select_file, index=False)


def fillin_err(
    select_files: list[os.PathLike], fillin_files: list[os.PathLike]
) -> None:
    metric_dfs = [
        pd.read_csv(select_file, header=0).assign(
            rpos1_d=lambda df: (df["rpos1_p"] - df["rpos1_e"]).abs(),
            rpos2_d=lambda df: (df["rpos2_p"] - df["rpos2_e"]).abs(),
            qpos1_d=lambda df: (df["qpos1_p"] - df["qpos1_e"]).abs(),
            qpos2_d=lambda df: (df["qpos2_p"] - df["qpos2_e"]).abs(),
        )
        for select_file in select_files
    ]
    worst_df = (
        pd
        .concat(metric_dfs, ignore_index=True)
        .groupby("idx", as_index=False)
        .agg(
            rpos1_d=pd.NamedAgg(column="rpos1_d", aggfunc="max"),
            rpos2_d=pd.NamedAgg(column="rpos2_d", aggfunc="max"),
            qpos1_d=pd.NamedAgg(column="qpos1_d", aggfunc="max"),
            qpos2_d=pd.NamedAgg(column="qpos2_d", aggfunc="max"),
        )
        .set_index("idx")
    )
    for metric_df, fillin_file in zip(metric_dfs, fillin_files):
        metric_df = metric_df[
            ["idx", "rpos1_d", "rpos2_d", "qpos1_d", "qpos2_d"]
        ].set_index("idx")
        fillin_df = (
            pd
            .concat(
                [
                    metric_df[col].combine_first(worst_df[col])
                    for col in ["rpos1_d", "rpos2_d", "qpos1_d", "qpos2_d"]
                ],
                axis=1,
            )
            .sort_index()
            .reset_index(names=["idx"])
        )
        pathlib.Path(os.fspath(fillin_file)).parent.mkdir(exist_ok=True, parents=True)
        fillin_df.to_csv(fillin_file, index=False)
