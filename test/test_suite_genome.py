#!/usr/bin/env python

import unittest

# Test cases depending on genome
from .test_get_sx_plasmid_file_ref.test_get_sx_plasmid_file_ref import (
    TestGetSxPlasmidFileRef,
)
from .test_work_flow.test_work_flow import TestWorkFlow


def suite_genome():
    suite = unittest.TestSuite()

    suite.addTest(TestGetSxPlasmidFileRef("test_get_sx_plasmid_file_ref"))
    suite.addTest(TestWorkFlow("test_work_flow"))

    return suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite_genome())
