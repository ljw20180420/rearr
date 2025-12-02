#!/usr/bin/env python

import unittest
import subprocess
from tempfile import mkstemp
import filecmp
import os


class TestGetSxPlasmidFileRef(unittest.TestCase):
    def setUp(self):
        self.plasmid_files = [
            "test/test_get_sx_plasmid_file_ref/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv",
            "test/test_get_sx_plasmid_file_ref/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G2.csv",
            "test/test_get_sx_plasmid_file_ref/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G3.csv",
            "test/test_get_sx_plasmid_file_ref/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A1.csv",
            "test/test_get_sx_plasmid_file_ref/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A2.csv",
            "test/test_get_sx_plasmid_file_ref/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A3.csv",
        ]

        self.ref_files = [
            mkstemp(dir="test/test_get_sx_plasmid_file_ref", suffix=".ref")[1]
            for _ in self.plasmid_files
        ]

    def test_get_sx_plasmid_file_ref(self):
        for plasmid_file, ref_file in zip(self.plasmid_files, self.ref_files):
            subprocess.run(
                f"""getSxPlasmidFileRef.sh {plasmid_file} {os.environ["GENOME"]} {os.environ["BOWTIE2INDEX"]} 50 0 10 100 > {ref_file}""",
                shell=True,
                executable="/bin/bash",
            )
            self.assertTrue(filecmp.cmp(ref_file, f"{plasmid_file}.ref", shallow=False))

    def tearDown(self):
        for ref_file in self.ref_files:
            os.remove(ref_file)


if __name__ == "__main__":
    unittest.main()
