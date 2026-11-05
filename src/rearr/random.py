from typing import Literal

import numpy as np
import pandas as pd


def random_DNA(length: int, rng: np.random.Generator) -> str:
    return "".join(rng.choice(["A", "C", "G", "T"], length))


def random_ref(length: int, cut: int, rng: np.random.Generator) -> str:
    DNA = rng.choice(["A", "C", "G", "T"], length)
    DNA[cut + 4] = "G"
    DNA[cut + 5] = "G"

    return "".join(DNA)


def SNP_DNA(DNA: str, error_rate: float, rng: np.random.Generator) -> str:
    DNA = np.array(list(DNA))
    SNP_mask = rng.random(len(DNA)) < error_rate
    DNA[SNP_mask] = rng.choice(["A", "C", "G", "T"], sum(SNP_mask))

    return "".join(DNA)


def indel_DNA(DNA: str, error_rate: float, rng: np.random.Generator) -> str:
    DNA = np.array(list(DNA))
    DNA[rng.random(len(DNA)) < error_rate] = ""
    inslens = rng.negative_binomial(1, 1 - error_rate, len(DNA) + 1)
    insDNAs = np.array([random_DNA(inslen, rng) for inslen in inslens])

    return "".join(np.stack([insDNAs, np.append(DNA, "")], axis=1).flatten())


def generate_ref(
    ref1len: int,
    ref2len: int,
    cut1: int,
    cut2: int,
    mode: Literal["single", "double"],
    ref_num: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    ref_df = pd.DataFrame({
        "ref1": [random_ref(ref1len, cut1, rng) for _ in range(ref_num)],
        "ref2": [random_ref(ref2len, cut2, rng) for _ in range(ref_num)],
    })

    if mode == "single":
        ref_df["ref1"] = ref_df["ref1"].str.slice(0, cut1) + ref_df["ref2"].str.slice(
            cut2, cut2 + ref1len - cut1
        )
        ref_df["ref2"] = ref_df["ref1"].str.slice(cut1 - cut2, cut1) + ref_df[
            "ref2"
        ].str.slice(cut2, None)

    return ref_df


def generate_edit(
    ref1: str,
    ref2: str,
    rstart1: int,
    rend1: int,
    rstart2: int,
    rend2: int,
    error_rate: float,
    query_per_ref: int,
    max_align_pos_shift_ratio: float,
    max_random_insertion_length: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    max_shift1 = int(len(ref1) * max_align_pos_shift_ratio)
    max_shift2 = int(len(ref2) * max_align_pos_shift_ratio)

    edit_df = pd.DataFrame({
        "rpos1": rend1 + rng.integers(-max_shift1, max_shift1 + 1, query_per_ref),
        "rpos2": rstart2 + rng.integers(-max_shift2, max_shift2 + 1, query_per_ref),
        "random_insertion": [
            random_DNA(random_insertion_length, rng)
            for random_insertion_length in rng.integers(
                low=0, high=max_random_insertion_length, size=query_per_ref
            )
        ],
    })

    edit_df = edit_df.assign(
        mut1=lambda df: df["rpos1"].map(
            lambda rpos1: indel_DNA(
                SNP_DNA(ref1[rstart1:rpos1], error_rate, rng), error_rate, rng
            )
        ),
        mut2=lambda df: df["rpos2"].map(
            lambda rpos2: indel_DNA(
                SNP_DNA(ref2[rpos2:rend2], error_rate, rng), error_rate, rng
            )
        ),
    )

    return edit_df
