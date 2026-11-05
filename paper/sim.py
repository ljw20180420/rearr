#!/usr/bin/env python

import os
import pathlib

from rearr import config, sim

if __name__ == "__main__":
    os.chdir(pathlib.Path(__file__).resolve().parent.parent)
    cfg = config.sim_config()

    sim.generate_ref(cfg)
    sim.generate_edit(cfg)

    sim.bench_rearr(cfg)
