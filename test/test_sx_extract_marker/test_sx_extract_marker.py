#!/usr/bin/env python

import unittest
import subprocess
from tempfile import mkstemp
import filecmp
import os


class TestSxExtractMarker(unittest.TestCase):
    def setUp(self):
        self.plasmid_files = [
            "test/test_sx_extract_marker/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv",
            "test/test_sx_extract_marker/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G2.csv",
            "test/test_sx_extract_marker/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G3.csv",
            "test/test_sx_extract_marker/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A1.csv",
            "test/test_sx_extract_marker/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A2.csv",
            "test/test_sx_extract_marker/final_hgsgrna_libb_all_0811_NAA_scaffold_nbt_A3.csv",
        ]

        self.marker1s = [
            mkstemp(dir="test/test_sx_extract_marker", suffix=".target.fa")[1]
            for _ in self.plasmid_files
        ]

        self.marker2s = [
            mkstemp(dir="test/test_sx_extract_marker", suffix=".pair.fa")[1]
            for _ in self.plasmid_files
        ]

    def test_sx_extract_marker(self):
        for plasmid_file, marker1, marker2 in zip(
            self.plasmid_files, self.marker1s, self.marker2s
        ):
            subprocess.run(
                f"""sxExtractMarker.sh {plasmid_file} > {marker1} 3> {marker2}""",
                shell=True,
                executable="/bin/bash",
            )
            self.assertTrue(
                filecmp.cmp(marker1, f"{plasmid_file}.target.fa", shallow=False)
            )
            self.assertTrue(
                filecmp.cmp(marker2, f"{plasmid_file}.pair.fa", shallow=False)
            )

    def tearDown(self):
        for marker1 in self.marker1s:
            os.remove(marker1)
        for marker2 in self.marker2s:
            os.remove(marker2)


if __name__ == "__main__":
    unittest.main()
