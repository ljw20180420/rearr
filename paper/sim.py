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
    # sim.bench_sim(cfg, lib="resso")
    # sim.bench_sim(cfg, lib="amplicon")
    # sim.bench_sim(cfg, lib="crisprvariants")
    # sim.bench_sim(cfg, lib="ampliconDIVider")
    # sim.bench_sim(cfg, lib="crispr_grant")
    # sim.bench_sim(cfg, lib="calculate_indel")
    # sim.bench_sim(cfg, lib="SelfTarget")
    # sim.bench_sim(cfg, lib="CRISPR_toolkit")
    # sim.bench_sim(cfg, lib="CRISPR-Analytics")
    # sim.bench_sim(cfg, lib="SIQ")
    sim.bench_sim(cfg, lib="PEM-Q")
