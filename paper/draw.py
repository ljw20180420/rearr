#!/usr/bin/env python

import os
import pathlib

from rearr import config, draw

if __name__ == "__main__":
    os.chdir(pathlib.Path(__file__).resolve().parent.parent)

    # draw.draw_sim_time_memory_benchmark(cfg=config.sim_config())
    # draw.draw_sim_accuracy_benchmark(cfg=config.sim_config())
    draw.draw_real_time_memory_benchmark(cfg=config.real_config())
