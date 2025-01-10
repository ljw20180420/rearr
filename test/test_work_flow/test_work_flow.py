#!/usr/bin/env python

import unittest
import subprocess
from tempfile import mkstemp
import filecmp
import os
import numpy as np
from ..utils.random_seq_methods import generate_random_DNA

class TestWorkFlow(unittest.TestCase):
    def setUp(self):
        self.toTestFiles = [
            *[
                f'test/test_work_flow/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv{ext}'
                for ext in ['.ref', '.correct']
            ],
            *[
                f'test/test_work_flow/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.{tp}.fa{ext}'
                for tp in ['target', 'pair']
                for ext in ["", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
            ],
            'test/test_work_flow/rearr.noDup',
            'test/test_work_flow/rearr.demultiplex',
            'test/test_work_flow/rearr.post',
            'test/test_work_flow/rearr.alg'
        ]

    def test_work_flow(self):
        subprocess.run(
            f'''./runWorkFlow.md -s''',
            shell=True,
            executable='/bin/bash'
        )
        for toTestFile in self.toTestFiles:
            self.assertTrue(
                filecmp.cmp(
                    toTestFile,
                    f'{toTestFile}.bak',
                    shallow=False
                )
            )
    
    def tearDown(self):
        for toTestFile in self.toTestFiles:
            os.remove(toTestFile)

if __name__ == '__main__':
    unittest.main()