#!/usr/bin/env python

import unittest
import subprocess
from tempfile import mkstemp
import filecmp
import os

class TestSxExtractSpliter(unittest.TestCase):
    def setUp(self):
        self.csvfiles = [
            'test/test_sx_extract_spliter/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv',
            'test/test_sx_extract_spliter/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G2.csv',
            'test/test_sx_extract_spliter/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G3.csv',
            'test/test_sx_extract_spliter/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A1.csv',
            'test/test_sx_extract_spliter/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A2.csv',
            'test/test_sx_extract_spliter/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A3.csv'
        ]

        self.spliter1s = [
            mkstemp(dir='test/test_sx_extract_spliter', suffix='.target.fa')[1]
            for _ in self.csvfiles
        ]

        self.spliter2s = [
            mkstemp(dir='test/test_sx_extract_spliter', suffix='.pair.fa')[1]
            for _ in self.csvfiles
        ]

    def test_sx_extract_spliter(self):
        for csvfile, spliter1, spliter2 in zip(self.csvfiles, self.spliter1s, self.spliter2s):
            subprocess.run(
                f'''sxExtractSpliter.md {csvfile} > {spliter1} 3> {spliter2}''',
                shell=True,
                executable='/bin/bash'
            )
            self.assertTrue(
                filecmp.cmp(
                    spliter1,
                    f'{csvfile}.target.fa',
                    shallow=False
                )
            )
            self.assertTrue(
                filecmp.cmp(
                    spliter2,
                    f'{csvfile}.pair.fa',
                    shallow=False
                )
            )

    def tearDown(self):
        for spliter1 in self.spliter1s:
            os.remove(spliter1)
        for spliter2 in self.spliter2s:
            os.remove(spliter2)


if __name__ == '__main__':
    unittest.main()