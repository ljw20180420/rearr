#!/usr/bin/env python

import os
import pathlib

from rearr import config, sim

if __name__ == "__main__":
    os.chdir(pathlib.Path(__file__).resolve().parent.parent)
    cfg = config.sim_config()

    # sim.generate_sim_ref(cfg)
    # sim.generate_sim_edit(cfg)
    # sim.bench_sim(cfg, lib="rearr")
    sim.bench_sim(cfg, lib="resso")
