#!/usr/bin/env python

import unittest

# Test cases not depending on genome
from .test_demultiplex.test_demultiplex import TestDemultiplex
from .test_rearr.test_rearr import TestRearrangementAlign, TestRearrangementHelp
from .test_remove_duplicates.test_remove_duplicates import TestRemoveDuplicates
from .test_work_flow.test_work_flow import TestWorkFlow


def suite_all():
    suite = unittest.TestSuite()

    suite.addTest(TestDemultiplex("test_demultiplex"))
    suite.addTest(TestRearrangementHelp("test_help"))
    suite.addTest(TestRearrangementAlign("test_align"))
    suite.addTest(TestRemoveDuplicates("test_remove_duplicates"))
    suite.addTest(TestWorkFlow("test_work_flow"))

    return suite


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite_all())
