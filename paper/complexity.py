#!/usr/bin/env python

import pandas as pd
import os
from plotnine import (
    ggplot,
    aes,
    geom_line,
    coord_equal,
    scale_x_continuous,
    scale_y_continuous,
)

os.chdir(os.path.dirname(os.path.abspath(__file__)))

os.makedirs("figures/complexity", exist_ok=True)


df = pd.read_csv("complexity.tsv", sep="\t")

for y in ["preprocess time", "total time", "align time", "DP time", "backtrack time"]:
    (
        ggplot(
            df,
            aes(
                "reference length",
                y,
                group="probability",
                color="probability",
            ),
        )
        + geom_line()
        + scale_x_continuous(trans="log2")
        + scale_y_continuous(trans="log2")
        + coord_equal()
    ).save(f"figures/complexity/{y}.png")
