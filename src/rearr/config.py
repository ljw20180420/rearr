import os
import pathlib


def complexity_config():
    return {
        "output_dir": pathlib.Path(os.environ["HOME"])
        / "sde1"
        / "rearr"
        / "complexity",
        "modes": ["single", "double"],
        "ref_num": 10,
        "ref_lens": [32, 64, 128, 256, 512, 1024, 2048, 4096],
        "error_rates": [0.01, 0.02, 0.03, 0.04, 0.05],
        "query_ratio": 0.5,
        "query_per_ref": 1000,
        "max_align_pos_shift_ratio": 0.05,
        "max_random_insertion_length": 10,
        "seed": 63036,
    }


def sim_config():
    pass
