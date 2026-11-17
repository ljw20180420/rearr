import pandas as pd
import plotly.graph_objects as go


def plotly_template() -> go.layout.Template:
    return go.layout.Template(
        layout=go.Layout(
            width=500,
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
    for target in ["wall_time", "memory"]:
        if target == "wall_time":
            draw_dir = cfg["output_dir"] / "draw" / "time"
        else:
            draw_dir = cfg["output_dir"] / "draw" / "memory"
        draw_dir.mkdir(exist_ok=True, parents=True)
        for mode in cfg["modes"]:
            fig = go.Figure(
                layout={
                    "width": 800,
                    "height": 500,
                }
            )
            for ref_len in cfg["ref_lens"]:
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

            fig.update_layout(
                scattermode="group",
                scattergap=0.5,
                legend_title_text="error rate",
                xaxis_title="package",
                yaxis_title="time (s)" if target == "wall_time" else "memory (kB)",
                template=plotly_template(),
            )
            fig.write_image(draw_dir / f"{mode}.pdf")


def draw_accuracy_benchmark(cfg) -> None:
    pass
