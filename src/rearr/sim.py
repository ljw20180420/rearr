import numpy as np
import pandas as pd

from .bench import get_bench_func
from .metric import fillin_err, select_best_indel
from .random import generate_edit, generate_ref


def generate_sim_ref(cfg: dict) -> None:
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
                ref_num=1,
                rng=rng,
            )

            (cfg["output_dir"] / "ref").mkdir(exist_ok=True, parents=True)
            ref_df.to_csv(
                cfg["output_dir"] / "ref" / f"{mode}_{ref_len}.csv",
                index=False,
            )


def generate_sim_edit(cfg: dict) -> None:
    rng = np.random.default_rng(cfg["seed"])

    for mode in cfg["modes"]:
        for ref_len in cfg["ref_lens"]:
            rstart = int(ref_len * (0.5 - cfg["query_ratio"] / 2))
            rend = int(ref_len * (0.5 + cfg["query_ratio"] / 2))

            ref_df = pd.read_csv(
                cfg["output_dir"] / "ref" / f"{mode}_{ref_len}.csv", header=0
            )

            for error_rate in cfg["error_rates"]:
                assert len(ref_df) == 1, "only single reference is allowed"
                edit_df = generate_edit(
                    ref1=ref_df.loc[0, "ref1"],
                    ref2=ref_df.loc[0, "ref2"],
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
                edit_df = edit_df.assign(ref_id=0)

                (cfg["output_dir"] / "edit").mkdir(exist_ok=True, parents=True)
                edit_df.to_csv(
                    cfg["output_dir"] / "edit" / f"{mode}_{ref_len}_{error_rate}.csv",
                    index=False,
                )


def bench_sim(cfg: dict, lib: str) -> None:
    (cfg["output_dir"] / lib / "predict").mkdir(exist_ok=True, parents=True)
    with open(cfg["output_dir"] / lib / "bench.csv", "w") as fd:
        fd.write("mode,ref_len,error_rate,wall_time,memory\n")
        fd.flush()
        for mode in cfg["modes"]:
            for ref_len in cfg["ref_lens"]:
                for error_rate in cfg["error_rates"]:
                    print(lib, mode, ref_len, error_rate)
                    wall_time, memory = get_bench_func(lib)(
                        ref_file=cfg["output_dir"] / "ref" / f"{mode}_{ref_len}.csv",
                        edit_file=cfg["output_dir"]
                        / "edit"
                        / f"{mode}_{ref_len}_{error_rate}.csv",
                        query_ratio=cfg["query_ratio"],
                        predict_file=cfg["output_dir"]
                        / lib
                        / "predict"
                        / f"{mode}_{ref_len}_{error_rate}.csv",
                    )

                    fd.write(f"{mode},{ref_len},{error_rate},{wall_time},{memory}\n")
                    fd.flush()


def metric_sim(cfg: dict, lib: str) -> None:
    (cfg["output_dir"] / lib / "metric").mkdir(exist_ok=True, parents=True)
    for mode in cfg["modes"]:
        for ref_len in cfg["ref_lens"]:
            for error_rate in cfg["error_rates"]:
                print(lib, mode, ref_len, error_rate)
                select_best_indel(
                    predict_file=cfg["output_dir"]
                    / lib
                    / "predict"
                    / f"{mode}_{ref_len}_{error_rate}.csv",
                    edit_file=cfg["output_dir"]
                    / "edit"
                    / f"{mode}_{ref_len}_{error_rate}.csv",
                    metric_file=cfg["output_dir"]
                    / lib
                    / "metric"
                    / f"{mode}_{ref_len}_{error_rate}.csv",
                )


def fillin_sim(cfg: dict) -> None:
    for mode in cfg["modes"]:
        for ref_len in cfg["ref_lens"]:
            for error_rate in cfg["error_rates"]:
                print(mode, ref_len, error_rate)
                fillin_err(
                    metric_files=[
                        cfg["output_dir"]
                        / lib
                        / "metric"
                        / f"{mode}_{ref_len}_{error_rate}.csv"
                        for lib in cfg["libs"]
                    ],
                    fillin_files=[
                        cfg["output_dir"]
                        / lib
                        / "fillin"
                        / f"{mode}_{ref_len}_{error_rate}.csv"
                        for lib in cfg["libs"]
                    ],
                )
