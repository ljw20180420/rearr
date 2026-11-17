import pandas as pd

from .bench import get_bench_func


def bench_real(cfg: dict, lib: str):
    (cfg["output_dir"] / lib / "predict").mkdir(exist_ok=True, parents=True)
    with open(cfg["output_dir"] / lib / "bench.csv", "w") as fd:
        fd.write("data_name,fq_file,wall_time,memory\n")
        fd.flush()
        for data_name, data_cfg in cfg["datasets"].items():
            fq2ref = pd.read_csv(
                data_cfg["fq2ref"], sep="\t", names=["fq_file", "ref_file"]
            ).sample(n=data_cfg["fq_num"], random_state=cfg["seed"], ignore_index=True)
            for fq_file, ref_file in fq2ref.itertuples(index=False):
                print(lib, data_name, fq_file, ref_file)
                ref_df = pd.read_csv(
                    data_cfg["ref_dir"] / ref_file,
                    sep="\t",
                    names=["rstart1", "ref1", "rend1", "rstart2", "ref2", "rend2"],
                )
                edit_df = (
                    pd
                    .DataFrame(
                        pd.read_csv(
                            data_cfg["fq_dir"] / f"{fq_file}.gz",
                            sep="\t",
                            header=None,
                            nrows=40 * data_cfg["sample_num"],
                        ).values.reshape(-1, 4),
                        columns=["seqname", "query", "plus", "quality"],
                    )
                    .assign(
                        query=lambda df: df["query"].str.upper(),
                    )[["query"]]
                    .drop_duplicates(ignore_index=True)
                    .head(n=data_cfg["sample_num"])
                    .assign(
                        ref_id=0,
                    )
                )
                wall_time, memory = get_bench_func(lib)(
                    ref_df=ref_df,
                    edit_df=edit_df,
                    predict_file=cfg["output_dir"] / lib / "predict" / f"{fq_file}.csv",
                )

                fd.write(f"{data_name},{fq_file},{wall_time},{memory}\n")
                fd.flush()
