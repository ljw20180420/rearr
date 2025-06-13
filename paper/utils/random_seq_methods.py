#!/usr/bin/env python
import numpy as np
import pandas as pd
import sys


class MutateSeq:
    def __init__(self, seed: int = 63036) -> None:
        self.rng: np.random.Generator = np.random.default_rng(seed)

    def random_DNA(self, length: int) -> str:
        return "".join(self.rng.choice(["A", "C", "G", "T"], length))

    def random_ref(self, length: int, cut: int) -> str:
        DNA = self.rng.choice(["A", "C", "G", "T"], length)
        DNA[cut + 4] = DNA[cut + 5] = "G"
        return "".join(DNA)

    def SNP_DNA(self, DNA: str, probability: float) -> str:
        DNA = np.array([DNA]).view("U1")
        SNP_mask = self.rng.random(len(DNA)) < probability
        DNA[SNP_mask] = self.rng.choice(["A", "C", "G", "T"], sum(SNP_mask))
        return "".join(DNA)

    def indel_DNA(self, DNA: str, probability: float) -> str:
        DNA = np.array([DNA]).view("U1")
        DNA[self.rng.random(len(DNA)) < probability] = ""
        inslens = self.rng.negative_binomial(1, 1 - probability, len(DNA) + 1)
        insDNAs = np.array([self.random_DNA(inslen) for inslen in inslens])
        return "".join(np.stack([insDNAs, np.append(DNA, "")]).ravel("F"))

    def CRISPR_editing(
        self,
        ref1: str,
        ref2: str,
        rstart1: int,
        rend1: int,
        rstart2: int,
        rend2: int,
        num: int,
        probability: float,
        shift_ratio: float = 0.05,
        rand_ins_range: list = [0, 10],
    ) -> pd.DataFrame:
        shift1 = int(len(ref1) * shift_ratio)
        shift2 = int(len(ref2) * shift_ratio)
        rpos1s = rend1 + self.rng.integers(-shift1, shift1 + 1, num)
        rpos2s = rstart2 + self.rng.integers(-shift2, shift2 + 1, num)
        inss = np.array(
            [
                self.random_DNA(rand_ins_len)
                for rand_ins_len in self.rng.integers(*rand_ins_range, num)
            ]
        )
        mut1s = np.array(
            [
                self.indel_DNA(
                    self.SNP_DNA(ref1[rstart1:rpos1], probability), probability
                )
                for rpos1 in rpos1s
            ]
        )
        mut2s = np.array(
            [
                self.indel_DNA(
                    self.SNP_DNA(ref2[rpos2:rend2], probability), probability
                )
                for rpos2 in rpos2s
            ]
        )
        qpos1s = np.char.str_len(mut1s)
        qpos2s = qpos1s + np.char.str_len(inss)

        return pd.DataFrame(
            {
                "query": np.char.add(
                    np.char.add(mut1s, inss),
                    mut2s,
                ),
                "count": [1] * num,
            }
        ), pd.DataFrame(
            {
                "rpos1": rpos1s,
                "rpos2": rpos2s,
                "qpos1": qpos1s,
                "qpos2": qpos2s,
            }
        )

    def rearr_API(
        self,
        ref_len: int,
        ref_num: int,
        query_per_ref: int,
        probability: float,
        ratio: float = 0.5,
        mode: str = "double",
    ):
        assert mode in ["single", "double"], "mode must be single or double"

        rstart = int(ref_len * (0.5 - ratio / 2))
        rend = int(ref_len * (0.5 + ratio / 2))

        ref_df = pd.DataFrame(
            {
                "rstart1": [rstart] * ref_num,
                "ref1": [self.random_ref(ref_len, rend) for _ in range(ref_num)],
                "rend1": [rend] * ref_num,
                "rstart2": [rstart] * ref_num,
                "ref2": [self.random_ref(ref_len, rstart) for _ in range(ref_num)],
                "rend2": [rend] * ref_num,
            }
        )
        if mode == "single":
            ref_df["ref1"] = ref_df["ref1"].str.slice(0, rend) + ref_df[
                "ref2"
            ].str.slice(rstart, rstart + ref_len - rend)
            ref_df["ref2"] = ref_df["ref1"].str.slice(rend - rstart, rend) + ref_df[
                "ref2"
            ].str.slice(rstart, None)
        query_dfs, truth_dfs = [], []
        for i in range(ref_num):
            query_df, truth_df = self.CRISPR_editing(
                ref1=ref_df.loc[i, "ref1"],
                ref2=ref_df.loc[i, "ref2"],
                rstart1=rstart,
                rend1=rend,
                rstart2=rstart,
                rend2=rend,
                num=query_per_ref,
                probability=probability,
            )
            query_df["ref_id"] = i
            query_dfs.append(query_df)
            truth_dfs.append(truth_df)

        return ref_df, pd.concat(query_dfs), pd.concat(truth_dfs)


if __name__ == "__main__":
    ref_len, ref_num, query_per_ref, probability, mode = sys.argv[1:]
    ref_len, ref_num, query_per_ref, probability = (
        int(ref_len),
        int(ref_num),
        int(query_per_ref),
        float(probability),
    )
    mutate_seq = MutateSeq()
    ref_df, query_df, truth_df = mutate_seq.rearr_API(
        ref_len, ref_num, query_per_ref, probability
    )
    ref_df.to_csv(sys.stdout, sep="\t", header=False, index=False)
    with open(3, "w") as fd:
        query_df.to_csv(fd, sep="\t", header=False, index=False)
    with open(4, "w") as fd:
        truth_df.to_csv(fd, sep="\t", header=False, index=False)
