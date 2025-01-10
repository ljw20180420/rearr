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

        # Write spliter seqs to .spliter files.
        self.spliter_num = 100
        self.spliters = [
            mkstemp(dir='test/test_demultiplex', suffix='.spliter')[1]
            for _ in range(3)
        ]
        for spliter in self.spliters:
            with open(spliter, 'w') as fd:
                for i in range(self.spliter_num):
                    spliter_seq = generate_random_DNA(rng.integers(50, 101), rng)
                    fd.write(f'>{i}\n{spliter_seq}\n')
            subprocess.run(
               f'''bowtie2-build -q {spliter} {spliter}''',
               shell=True,
               executable='/bin/bash'
            )

        # Generate .noDup file.
        self.noDupNum = 1000
        self.rmDupFile = mkstemp(dir='test/test_demultiplex', suffix='.noDup')[1]
        with open(self.rmDupFile, 'w') as wd:
            spliters_seqs = []
            for spliter in self.spliters:
                with open(spliter, 'r') as fd:
                    spliters_seqs.append(
                        fd.read().splitlines()[1::2]
                    )
            for i in range(self.noDupNum):
                ridx = rng.integers(self.spliter_num)
                r_spliter_seqs = [
                    spliter_seqs[ridx]
                    for spliter_seqs in spliters_seqs
                ]
                for r_spliter_seq in r_spliter_seqs:
                    if rng.random() < 0.1:
                        noDup_seq = generate_random_DNA(rng.integers(50, 101), rng)
                    else:
                        s_start = rng.integers(11)
                        s_end = len(r_spliter_seq) - rng.integers(11)
                        noDup_seq = r_spliter_seq[s_start:s_end]
                        noDup_seq = SNP_DNA(noDup_seq, 0.01, rng)
                        noDup_seq = indel_DNA(noDup_seq, 0.01, rng)
                        noDup_seq = ''.join(
                            [
                                generate_random_DNA(rng.integers(11), rng),
                                noDup_seq,
                                generate_random_DNA(rng.integers(11), rng)
                            ]
                        )
                    wd.write(f'{noDup_seq}\t')
                wd.write(f'{rng.integers(11)}\n')
                    
        # Set minimal scores.
        self.minScores = [str(30) for _ in self.spliters]

        self.demultiplexFile = mkstemp(dir='test/test_demultiplex', suffix='.demultiplex')[1]

    def test_demultiplex(self):
        subprocess.run(
            f'''spliterIndices={','.join(self.spliters)} minScores={','.join(self.minScores)} demultiplex.md {self.rmDupFile} > {self.demultiplexFile}''',
            shell=True,
            executable='/bin/bash'
        )
        self.assertTrue(
            filecmp.cmp(
                self.demultiplexFile,
                'test/test_demultiplex/output.demultiplex',
                shallow=False
            )
        )

    def tearDown(self):
        for spliter in self.spliters:
            for ext in ['', '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']:
                os.remove(f'{spliter}{ext}')
        os.remove(self.rmDupFile)
        os.remove(self.demultiplexFile)

if __name__ == '__main__':
    unittest.main()

# spliterIndices=index1,index2,... minScores=score1,score2,... demultiplex.md inputFile >demultiplexFile