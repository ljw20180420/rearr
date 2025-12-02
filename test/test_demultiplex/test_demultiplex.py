#!/usr/bin/env python

import unittest
import subprocess
import numpy as np
from tempfile import mkstemp
import filecmp
import os
from ..utils.random_seq_methods import generate_random_DNA, SNP_DNA, indel_DNA


class TestDemultiplex(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(63036)

        # Write marker seqs to .marker files.
        self.marker_num = 100
        self.markers = [
            mkstemp(dir="test/test_demultiplex", suffix=".marker")[1] for _ in range(3)
        ]
        for marker in self.markers:
            with open(marker, "w") as fd:
                for i in range(self.marker_num):
                    marker_seq = generate_random_DNA(rng.integers(50, 101), rng)
                    fd.write(f">{i}\n{marker_seq}\n")
            subprocess.run(
                f"""bowtie2-build -q {marker} {marker}""",
                shell=True,
                executable="/bin/bash",
            )

        # Generate .noDup file.
        self.noDupNum = 1000
        self.rmDupFile = mkstemp(dir="test/test_demultiplex", suffix=".noDup")[1]
        with open(self.rmDupFile, "w") as wd:
            markers_seqs = []
            for marker in self.markers:
                with open(marker, "r") as fd:
                    markers_seqs.append(fd.read().splitlines()[1::2])
            for i in range(self.noDupNum):
                ridx = rng.integers(self.marker_num)
                r_marker_seqs = [marker_seqs[ridx] for marker_seqs in markers_seqs]
                for r_marker_seq in r_marker_seqs:
                    if rng.random() < 0.1:
                        noDup_seq = generate_random_DNA(rng.integers(50, 101), rng)
                    else:
                        s_start = rng.integers(11)
                        s_end = len(r_marker_seq) - rng.integers(11)
                        noDup_seq = r_marker_seq[s_start:s_end]
                        noDup_seq = SNP_DNA(noDup_seq, 0.01, rng)
                        noDup_seq = indel_DNA(noDup_seq, 0.01, rng)
                        noDup_seq = "".join(
                            [
                                generate_random_DNA(rng.integers(11), rng),
                                noDup_seq,
                                generate_random_DNA(rng.integers(11), rng),
                            ]
                        )
                    wd.write(f"{noDup_seq}\t")
                wd.write(f"{rng.integers(11)}\n")

        # Set minimal scores.
        self.minScores = [str(30) for _ in self.markers]

        self.demultiplexFile = mkstemp(
            dir="test/test_demultiplex", suffix=".demultiplex"
        )[1]

    def test_demultiplex(self):
        subprocess.run(
            f"""markerIndices={','.join(self.markers)} minScores={','.join(self.minScores)} demultiplex.sh {self.rmDupFile} > {self.demultiplexFile}""",
            shell=True,
            executable="/bin/bash",
        )
        self.assertTrue(
            filecmp.cmp(
                self.demultiplexFile,
                "test/test_demultiplex/output.demultiplex",
                shallow=False,
            )
        )

    def tearDown(self):
        for marker in self.markers:
            for ext in [
                "",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ]:
                os.remove(f"{marker}{ext}")
        os.remove(self.rmDupFile)
        os.remove(self.demultiplexFile)


if __name__ == "__main__":
    unittest.main()
