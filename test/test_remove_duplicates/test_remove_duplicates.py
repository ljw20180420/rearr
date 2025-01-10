#!/usr/bin/env python

import unittest
import subprocess
import numpy as np
from tempfile import mkstemp
import filecmp
import os
from ..utils.random_seq_methods import generate_random_DNA

class TestRemoveDuplicates(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(63036)

        self.fastqs = [
            mkstemp(dir='test/test_remove_duplicates', suffix=rng.choice(['.fq', '.fastq']))[1]
            for _ in range(3)
        ]
        self.seq_num = 1000
        self.rmDupFile = mkstemp(dir='test/test_remove_duplicates', suffix='.noDup')[1]

        idxs = rng.permuted(
            np.arange(self.seq_num).repeat(rng.integers(1, 11, self.seq_num))
        )
        for fastq in self.fastqs:
            with open(fastq, 'w') as fd:
                seqs = [
                    generate_random_DNA(rng.integers(50, 101), rng)
                    for _ in range(self.seq_num)
                ]       
                for i, idx in enumerate(idxs):
                    seq = seqs[idx]
                    qual = "".join(
                        [
                            chr(q33)
                            for q33 in rng.integers(33, 126, len(seq))
                        ]
                    )
                    fd.write(f'''@seq{i}\n{seq}\n+\n{qual}\n''')

    def test_remove_duplicates(self):
        subprocess.run(
            f'''removeDuplicates.md {' '.join(self.fastqs)} > {self.rmDupFile}''',
            shell=True,
            executable='/bin/bash'
        )
        self.assertTrue(
            filecmp.cmp(
                self.rmDupFile,
                'test/test_remove_duplicates/output.noDup',
                shallow=False
            )
        )

    def tearDown(self):
        for fastq in self.fastqs:
            os.remove(fastq)
        os.remove(self.rmDupFile)

if __name__ == '__main__':
    unittest.main()