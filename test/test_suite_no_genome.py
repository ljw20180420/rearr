#!/usr/bin/env python

import unittest
# Test cases not depending on genome
from .test_demultiplex.test_demultiplex import TestDemultiplex
from .test_rearr.test_rearr import TestRearrangementHelp, TestRearrangementAlign
from .test_remove_duplicates.test_remove_duplicates import TestRemoveDuplicates
from .test_sx_cut_r2_adapter_filter_cumulate.test_sx_cut_r2_adapter_filter_cumulate import TestSxCumulateToMapCutAdaptSpliter, TestSxCutR2AdapterFilterCumulate
from .test_sx_extract_spliter.test_sx_extract_spliter import TestSxExtractSpliter

def suite_no_genome():
    suite = unittest.TestSuite()

    suite.addTest(TestDemultiplex('test_demultiplex'))
    suite.addTest(TestRearrangementHelp('test_help'))
    suite.addTest(TestRearrangementAlign('test_align'))
    suite.addTest(TestRemoveDuplicates('test_remove_duplicates'))
    suite.addTest(TestSxCumulateToMapCutAdaptSpliter('test_sx_cumulate_to_map_cut_adapt_spliter'))
    suite.addTest(TestSxCutR2AdapterFilterCumulate('test_sx_cut_r2_adapter_filter_cumulate'))
    suite.addTest(TestSxExtractSpliter('test_sx_extract_spliter'))

    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite_no_genome())