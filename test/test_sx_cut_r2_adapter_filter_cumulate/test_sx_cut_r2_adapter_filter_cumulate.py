#!/usr/bin/env python

import unittest
import subprocess
from tempfile import mkstemp
import filecmp
import os
import numpy as np
from ..utils.random_seq_methods import generate_random_DNA

class TestSxCumulateToMapCutAdaptSpliter(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(63036)

        # Generate file with consecutive duplicates.
        self.toAccFile = mkstemp(dir='test/test_sx_cut_r2_adapter_filter_cumulate', suffix='.toAcc')[1]
        with open(self.toAccFile, 'w') as fd:
            for _ in range(1000):
                seq = generate_random_DNA(rng.integers(50, 101), rng)
                ref_id = rng.integers(100)
                for repe in range(rng.integers(11)):
                    fd.write(f'{seq}\t{rng.integers(11)}\t{ref_id}\n')

        self.postFilePartialTest = mkstemp(dir='test/test_sx_cut_r2_adapter_filter_cumulate', suffix='.post.partial.test')[1]

    def test_sx_cumulate_to_map_cut_adapt_spliter(self):
        subprocess.run(
            f'''gawk -f sxCumulateToMapCutAdaptSpliter.awk < {self.toAccFile} > {self.postFilePartialTest}''',
            shell=True,
            executable='/bin/bash'
        )
        self.assertTrue(
            filecmp.cmp(
                self.postFilePartialTest,
                'test/test_sx_cut_r2_adapter_filter_cumulate/output.post.partial.test',
                shallow=False
            )
        )
                
    def tearDown(self):
        os.remove(self.postFilePartialTest)
        os.remove(self.toAccFile)

class TestSxCutR2AdapterFilterCumulate(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(63036)

        self.demultiplexFile = mkstemp(dir='test/test_sx_cut_r2_adapter_filter_cumulate', suffix='.demultiplex')[1]
        with open(self.demultiplexFile, 'w') as fd:
            for _ in range(1000):
                seq1 = generate_random_DNA(rng.integers(50, 101), rng)
                seq2 = generate_random_DNA(rng.integers(50, 101), rng)
                ref_id = rng.integers(100)
                # spliter_start1, spliter_end1, seq_start1, seq_end1, spliter_start2, spliter_end2, seq_start2, seq_end2
                spliter_poses = rng.integers(0, 11, 8)
                for repe in range(rng.integers(11)):
                    fd.write(f'{seq1}\t{seq2}\t{rng.integers(11)}\t{ref_id}')
                    for spliter_pos in spliter_poses:
                        fd.write(f'\t{spliter_pos}')
                    fd.write('\n')

        self.postFile = mkstemp(dir='test/test_sx_cut_r2_adapter_filter_cumulate', suffix='.post')[1]

    def test_sx_cut_r2_adapter_filter_cumulate(self):
        subprocess.run(
            f'''sxCutR2AdapterFilterCumulate.md {self.demultiplexFile} 30 > {self.postFile}''',
            shell=True,
            executable='/bin/bash'
        )
        self.assertTrue(
            filecmp.cmp(
                self.postFile,
                'test/test_sx_cut_r2_adapter_filter_cumulate/output.post',
                shallow=False
            )
        )

    def tearDown(self):
        os.remove(self.demultiplexFile)
        os.remove(self.postFile)

if __name__ == '__main__':
    unittest.main()