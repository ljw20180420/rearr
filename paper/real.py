#!/usr/bin/env python

import os
import pathlib

from rearr import config, draw, real

if __name__ == "__main__":
    os.chdir(pathlib.Path(__file__).resolve().parent.parent)
    cfg = config.real_config()

    for lib in cfg["libs"]:
        real.bench_real(cfg, lib=lib)

    draw.draw_real_time_memory_benchmark(cfg=cfg)
