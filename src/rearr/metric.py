import os

import pandas as pd


def select_best_indel(
    predict_file: os.PathLike, edit_file: os.PathLike, metric_file: os.PathLike
) -> pd.DataFrame:
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
