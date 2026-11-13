import os
import pathlib

import pandas as pd


def select_best_indel(
    predict_file: os.PathLike, edit_file: os.PathLike, metric_file: os.PathLike
) -> None:
    predict_df = pd.read_csv(predict_file, header=0)
    edit_df = pd.read_csv(edit_file, header=0, keep_default_na=False)

    edit_df = edit_df.assign(
        qpos1=lambda df: df["mut1"].str.len(),
        qpos2=lambda df: df["qpos1"] + df["random_insertion"].str.len(),
    ).reset_index(names=["idx"])

    predict_df = predict_df.merge(
        right=edit_df,
        how="left",
        on="idx",
        suffixes=("_p", "_e"),
        validate="many_to_one",
    )
    predict_df = predict_df.assign(
        total_err=lambda df: (
            (df["rpos1_p"] - df["rpos1_e"]).abs()
            + (df["rpos2_p"] - df["rpos2_e"]).abs()
            + (df["qpos1_p"] - df["qpos1_e"]).abs()
            + (df["qpos2_p"] - df["qpos2_e"]).abs()
        ),
        min_total_err=lambda df: df.groupby("idx")["total_err"].transform("min"),
    )
    predict_df = predict_df.query("total_err == min_total_err").reset_index(drop=True)
    predict_df = predict_df.assign(
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
    predict_df = predict_df.query("center_err == min_center_err").reset_index(drop=True)
    predict_df = predict_df.groupby("idx", as_index=False).first()
    predict_df[
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
    ].to_csv(metric_file, index=False)


def fillin_err(
    metric_files: list[os.PathLike], fillin_files: list[os.PathLike]
) -> None:
    metric_dfs = [
        pd.read_csv(metric_file, header=0).assign(
            rpos1_d=lambda df: (df["rpos1_p"] - df["rpos1_e"]).abs(),
            rpos2_d=lambda df: (df["rpos2_p"] - df["rpos2_e"]).abs(),
            qpos1_d=lambda df: (df["qpos1_p"] - df["qpos1_e"]).abs(),
            qpos2_d=lambda df: (df["qpos2_p"] - df["qpos2_e"]).abs(),
        )
        for metric_file in metric_files
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
