import numpy as np
import pandas as pd
import plotly.express as px

from .bench import bench_rearr_full
from .random import generate_edit, generate_ref


def generate_complexity_ref(cfg: dict) -> None:
    rng = np.random.default_rng(cfg["seed"])

    for mode in cfg["modes"]:
        for ref_len in cfg["ref_lens"]:
            rstart = int(ref_len * (0.5 - cfg["query_ratio"] / 2))
            rend = int(ref_len * (0.5 + cfg["query_ratio"] / 2))

            ref_df = generate_ref(
                ref1len=ref_len,
                ref2len=ref_len,
                cut1=rend,
                cut2=rstart,
                mode=mode,
                ref_num=cfg["ref_num"],
                rng=rng,
            )

            (cfg["output_dir"] / "ref").mkdir(exist_ok=True, parents=True)
            ref_df.to_csv(
                cfg["output_dir"] / "ref" / f"{mode}_{ref_len}.csv",
                index=False,
            )


def generate_complexity_edit(cfg: dict) -> None:
    rng = np.random.default_rng(cfg["seed"])

    for mode in cfg["modes"]:
        for ref_len in cfg["ref_lens"]:
            rstart = int(ref_len * (0.5 - cfg["query_ratio"] / 2))
            rend = int(ref_len * (0.5 + cfg["query_ratio"] / 2))

            ref_df = pd.read_csv(
                cfg["output_dir"] / "ref" / f"{mode}_{ref_len}.csv", header=0
            )

            for error_rate in cfg["error_rates"]:
                edit_dfs = []
                for i, ref1, ref2 in ref_df[["ref1", "ref2"]].itertuples():
                    edit_df = generate_edit(
                        ref1=ref1,
                        ref2=ref2,
                        rstart1=rstart,
                        rend1=rend,
                        rstart2=rstart,
                        rend2=rend,
                        error_rate=error_rate,
                        query_per_ref=cfg["query_per_ref"],
                        max_align_pos_shift_ratio=cfg["max_align_pos_shift_ratio"],
                        max_random_insertion_length=cfg["max_random_insertion_length"],
                        rng=rng,
                    )
                    edit_df = edit_df.assign(ref_id=i)
                    edit_dfs.append(edit_df)

                (cfg["output_dir"] / "edit").mkdir(exist_ok=True, parents=True)
                pd.concat(edit_dfs, ignore_index=True).to_csv(
                    cfg["output_dir"] / "edit" / f"{mode}_{ref_len}_{error_rate}.csv",
                    index=False,
                )


def bench_complexity_rearr(cfg: dict) -> None:
    (cfg["output_dir"] / "predict").mkdir(exist_ok=True, parents=True)
    with open(cfg["output_dir"] / "bench.csv", "w") as fd:
        fd.write(
            "mode,ref_len,error_rate,total_time,preprocess_time,align_time,DP_time,backtrack_time,wall_time,memory\n"
        )
        fd.flush()
        for mode in cfg["modes"]:
            for ref_len in cfg["ref_lens"]:
                for error_rate in cfg["error_rates"]:
                    print(mode, ref_len, error_rate)
                    (
                        total_time,
                        preprocess_time,
                        align_time,
                        DP_time,
                        backtrack_time,
                        wall_time,
                        memory,
                    ) = bench_rearr_full(
                        ref_file=cfg["output_dir"] / "ref" / f"{mode}_{ref_len}.csv",
                        edit_file=cfg["output_dir"]
                        / "edit"
                        / f"{mode}_{ref_len}_{error_rate}.csv",
                        query_ratio=cfg["query_ratio"],
                        predict_file=cfg["output_dir"]
                        / "predict"
                        / f"{mode}_{ref_len}_{error_rate}.predict",
                    )
                    fd.write(
                        f"{mode},{ref_len},{error_rate},{total_time},{preprocess_time},{align_time},{DP_time},{backtrack_time},{wall_time},{memory}\n"
                    )
                    fd.flush()


def draw_complexity_rearr(cfg: dict) -> None:
    (cfg["output_dir"] / "draw").mkdir(exist_ok=True, parents=True)
    df = pd.read_csv(cfg["output_dir"] / "bench.csv", header=0)
    for mode in cfg["modes"]:
        df_slice = df.query("mode == @mode").reset_index(drop=True)
        for y in [
            "total_time",
            "preprocess_time",
            "align_time",
            "DP_time",
            "backtrack_time",
        ]:
            fig = px.line(
                df_slice,
                x="ref_len",
                y=y,
                color="error_rate",
                log_x=True,
                log_y=True,
            )
            fig.update_yaxes(scaleanchor="x", scaleratio=1)
            fig.update_traces(line={"width": 0.5})
            fig.write_image(cfg["output_dir"] / "draw" / f"{mode}_{y}.pdf")
