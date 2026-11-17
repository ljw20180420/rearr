import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def plotly_template() -> go.layout.Template:
    return go.layout.Template(
        layout=go.Layout(
            width=800,
            height=500,
            font={
                "size": 12,
                "color": "black",
            },
            xaxis={
                "ticks": "outside",
                "ticklen": 5,
                "automargin": True,
            },
            yaxis={
                "ticks": "outside",
                "ticklen": 5,
                "automargin": True,
            },
        )
    )


def draw_time_memory_benchmark(cfg) -> None:
    color_dict = {
        0.01: "red",
        0.02: "green",
        0.03: "blue",
        0.04: "purple",
        0.05: "black",
    }
    df = pd.concat(
        [
            pd.read_csv(cfg["output_dir"] / lib / "bench.csv", header=0).assign(lib=lib)
            for lib in cfg["libs"]
        ],
        ignore_index=True,
    )
    scattergap = 0.3
    offsets = [
        (scattergap - 1) / 2 + (ri + 0.5) / len(cfg["ref_lens"]) * (1 - scattergap)
        for ri in range(len(cfg["ref_lens"]))
    ]
    for target in ["wall_time", "memory"]:
        if target == "wall_time":
            draw_dir = cfg["output_dir"] / "draw" / "time"
        else:
            draw_dir = cfg["output_dir"] / "draw" / "memory"
        draw_dir.mkdir(exist_ok=True, parents=True)
        for mode in cfg["modes"]:
            fig = go.Figure()
            for ri, ref_len in enumerate(cfg["ref_lens"]):
                for error_rate in cfg["error_rates"]:
                    df_slice = df.query(
                        "mode == @mode and ref_len == @ref_len and error_rate == @error_rate"
                    ).reset_index(drop=True)

                    fig = fig.add_trace(
                        go.Scatter(
                            x=df_slice["lib"],
                            y=df_slice[target],
                            offsetgroup=ref_len,
                            mode="markers",
                            marker={
                                "symbol": "line-ew",
                                "size": 2,
                                "line": {
                                    "width": 1,
                                    "color": color_dict[error_rate],
                                },
                            },
                            name=error_rate,
                            legendgroup=error_rate,
                            showlegend=ref_len == 100,
                        )
                    )
                for x in range(len(cfg["libs"])):
                    fig.add_annotation(
                        x=x + offsets[ri],
                        y=0,
                        text=str(ref_len),
                        showarrow=False,
                        textangle=90,
                        xanchor="center",
                        yanchor="top",
                        font={
                            "size": 8,
                        },
                    )

            fig.update_xaxes(
                type="category", categoryorder="array", categoryarray=cfg["libs"]
            )
            fig.update_layout(
                scattermode="group",
                scattergap=scattergap,
                legend_title_text="error rate",
                xaxis_title="package",
                yaxis_title="time (s)" if target == "wall_time" else "memory (kB)",
                template=plotly_template(),
            )
            fig.write_image(draw_dir / f"{mode}.pdf")


def draw_accuracy_benchmark(cfg) -> None:
    error_name_dict = {
        "rpos1_d": "reference upstream",
        "rpos2_d": "reference downstream",
        "qpos1_d": "query upstream",
        "qpos2_d": "query downstream",
    }
    color_dict = {
        "rpos1_d": "red",
        "rpos2_d": "green",
        "qpos1_d": "blue",
        "qpos2_d": "purple",
    }
    (cfg["output_dir"] / "draw" / "accuracy").mkdir(exist_ok=True, parents=True)
    for mode in cfg["modes"]:
        for ref_len in cfg["ref_lens"]:
            for error_rate in cfg["error_rates"]:
                df = (
                    pd
                    .concat(
                        [
                            pd.read_csv(
                                cfg["output_dir"]
                                / lib
                                / "fillin"
                                / f"{mode}_{ref_len}_{error_rate}.csv",
                                header=0,
                            ).assign(lib=lib)
                            for lib in cfg["libs"]
                        ],
                        ignore_index=True,
                    )
                    .groupby("lib", as_index=False)
                    .agg(
                        rpos1_d_mean=pd.NamedAgg(
                            column="rpos1_d",
                            aggfunc="mean",
                        ),
                        rpos1_d_sem=pd.NamedAgg(
                            column="rpos1_d",
                            aggfunc="sem",
                        ),
                        rpos2_d_mean=pd.NamedAgg(
                            column="rpos2_d",
                            aggfunc="mean",
                        ),
                        rpos2_d_sem=pd.NamedAgg(
                            column="rpos2_d",
                            aggfunc="sem",
                        ),
                        qpos1_d_mean=pd.NamedAgg(
                            column="qpos1_d",
                            aggfunc="mean",
                        ),
                        qpos1_d_sem=pd.NamedAgg(
                            column="qpos1_d",
                            aggfunc="sem",
                        ),
                        qpos2_d_mean=pd.NamedAgg(
                            column="qpos2_d",
                            aggfunc="mean",
                        ),
                        qpos2_d_sem=pd.NamedAgg(
                            column="qpos2_d",
                            aggfunc="sem",
                        ),
                    )
                )
                fig = go.Figure()
                for target in ["rpos1_d", "rpos2_d", "qpos1_d", "qpos2_d"]:
                    fig = fig.add_trace(
                        go.Scatter(
                            x=df["lib"],
                            y=df[f"{target}_mean"],
                            offsetgroup=target,
                            mode="markers",
                            marker={
                                "color": color_dict[target],
                            },
                            name=error_name_dict[target],
                            error_y={
                                "type": "data",
                                "array": df[f"{target}_sem"],
                                "visible": True,
                            },
                        ),
                    )

                fig.update_xaxes(
                    type="category", categoryorder="array", categoryarray=cfg["libs"]
                )
                fig.update_layout(
                    scattermode="group",
                    scattergap=0.4,
                    legend_title_text="junction error type",
                    xaxis_title="package",
                    yaxis_title="error (bp)",
                    template=plotly_template(),
                )
                fig.write_image(
                    cfg["output_dir"]
                    / "draw"
                    / "accuracy"
                    / f"{mode}_{ref_len}_{error_rate}.pdf"
                )
