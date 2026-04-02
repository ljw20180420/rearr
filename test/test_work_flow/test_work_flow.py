#!/usr/bin/env python

import filecmp
import os
import subprocess
import unittest


class TestWorkFlow(unittest.TestCase):
    def setUp(self):
        self.toTestFiles = [
            *[
                f"{tp}.fa.{suffix}"
                for tp in ["target", "pair"]
                for suffix in [
                    "1.bt2",
                    "2.bt2",
                    "3.bt2",
                    "4.bt2",
                    "rev.1.bt2",
                    "rev.2.bt2",
                ]
            ],
            "test/test_work_flow/ref.direct",
            "test/test_work_flow/rearr.noDup",
            "test/test_work_flow/rearr.demultiplex",
            "test/test_work_flow/rearr.post",
            "test/test_work_flow/rearr.alg",
        ]

    def test_work_flow(self):
        subprocess.run(f"""./runWorkFlow.sh -s""", shell=True, executable="/bin/bash")

        for toTestFile in self.toTestFiles:
            self.assertTrue(filecmp.cmp(toTestFile, f"{toTestFile}.bak", shallow=False))

    def tearDown(self):
        for toTestFile in self.toTestFiles:
            os.remove(toTestFile)


if __name__ == "__main__":
    unittest.main()
