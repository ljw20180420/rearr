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
    return {
        "output_dir": pathlib.Path(os.environ["HOME"]) / "sde1" / "rearr" / "sim",
        "modes": ["single", "double"],
        "ref_lens": [100, 200, 300, 400, 500],
        "error_rates": [0.01, 0.02, 0.03, 0.04, 0.05],
        "query_ratio": 0.5,
        "query_per_ref": 100,
        "max_align_pos_shift_ratio": 0.05,
        "max_random_insertion_length": 10,
        "seed": 63036,
        "libs": [
            "rearr",
            "resso",
            "amplicon",
            "crisprvariants",
            "ampliconDIVider",
            "crispr_grant",
            "calculate_indel",
            "SelfTarget",
            "CRISPR_toolkit",
            "CRISPR-Analytics",
            "SIQ",
            "PEM-Q",
        ],
    }


def real_config():
    return {
        "output_dir": pathlib.Path(os.environ["HOME"]) / "sde1" / "rearr" / "real",
        "seed": 63036,
        "datasets": {
            "SJLJH": {
                "fq_num": 10,
                "sample_num": 1000,
                "fq_dir": pathlib.Path(os.environ["HOME"])
                / "sdc1"
                / "sx_data"
                / "SJLJH"
                / "data",
                "fq2ref": pathlib.Path(os.environ["HOME"])
                / "new_fold"
                / "projects"
                / "CRISPR_align"
                / "SJLJH"
                / "fq2ref.tsv",
                "ref_dir": pathlib.Path(os.environ["HOME"])
                / "sdc1"
                / "sx_data"
                / "SJLJH"
                / "refs",
            },
        },
        "libs": [
            "rearr",
            "resso",
            "amplicon",
            "crisprvariants",
            "ampliconDIVider",
            "crispr_grant",
            "calculate_indel",
            "SelfTarget",
            "CRISPR_toolkit",
            "CRISPR-Analytics",
            "SIQ",
            "PEM-Q",
        ],
    }
