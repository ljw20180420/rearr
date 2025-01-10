#!/usr/bin/env python

import unittest
import subprocess
import numpy as np
from tempfile import mkstemp
import filecmp
import os
from ..utils.random_seq_methods import generate_random_DNA, SNP_DNA, indel_DNA

class TestRearrangementHelp(unittest.TestCase):
    def test_help(self):
        result = subprocess.run("rearrangement -h", capture_output=True, shell=True, executable="/bin/bash")
        self.assertTrue(result.stdout.decode().startswith('###Basic Usage'))
        result = subprocess.run("rearrangement -help", capture_output=True, shell=True, executable="/bin/bash")
        self.assertTrue(result.stdout.decode().startswith('###Basic Usage'))
        result = subprocess.run("rearrangement --help", capture_output=True, shell=True, executable="/bin/bash")
        self.assertTrue(result.stdout.decode().startswith('###Basic Usage'))

class TestRearrangementAlign(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(63036)

        # Write random references to temp .ref file. Write random correct directions to temp .correct file.
        _, self.ref_file = mkstemp(dir="test/test_rearr", suffix=".ref")
        _, self.correct_file = mkstemp(dir="test/test_rearr", suffix=".correct")
        with open(self.ref_file, 'w') as ref_fd, open(self.correct_file, 'w') as correct_fd:
            for seg_num in rng.integers(1, 5, 100):
                for seg_idx in range(seg_num):
                    ref_len = rng.integers(50, 101)
                    ref = generate_random_DNA(ref_len, rng)
                    cut1 = rng.integers(0, 11)
                    cut2 = ref_len - rng.integers(0, 11)
                    sep = '\t' if seg_idx < seg_num - 1 else '\n'
                    _ = ref_fd.write(f'{cut1}\t{ref}\t{cut2}{sep}')
                    if seg_idx < seg_num - 1:
                        correct_type = 'up' if rng.random() < 0.5 else 'down'
                        correct_sep = '\t' if seg_idx < seg_num - 2 else ''
                        _ = correct_fd.write(f'{correct_type}{correct_sep}')
                correct_fd.write('\n')

        # Read random references into ref_lines.
        with open(self.ref_file, 'r') as fd:
            ref_lines = []
            for line in fd:
                fields = line.strip().split()
                ref_line = {
                    'cut1': [],
                    'cut2': [],
                    'ref': []
                }
                for i in range(0, len(fields), 3):
                    ref_line['cut1'].append(int(fields[i]))
                    ref_line['cut2'].append(int(fields[i+2]))
                    ref_line['ref'].append(fields[i+1])
                ref_lines.append(ref_line)

        # Write random query sequences to temp .post file.
        self.query_num = 1000
        _, self.post_file = mkstemp(dir="test/test_rearr", suffix=".post")
        with open(self.post_file, 'w') as fd:
            for ref_id in rng.integers(0, len(ref_lines), self.query_num):
                ref_line = ref_lines[ref_id]
                refs, cut1s, cut2s = ref_line['ref'], ref_line['cut1'], ref_line['cut2']
                query_segs = []
                for ref, cut1, cut2 in zip(refs, cut1s, cut2s):
                    query_segs.append(generate_random_DNA(rng.integers(4), rng))
                    ref_start = rng.integers(cut1 + 11)
                    ref_end = len(ref) - rng.integers(len(ref) - cut2 + 11)
                    query_seg = ref[ref_start:ref_end]
                    query_seg = SNP_DNA(query_seg, 0.02, rng)
                    query_seg = indel_DNA(query_seg, 0.02, rng)
                    query_segs.append(query_seg)
                    query_segs.append(generate_random_DNA(rng.integers(4), rng))
                _ = fd.write(f'''{"".join(query_segs)}\t{rng.integers(100)}\t{ref_id}\n''')

        _, self.alg_file = mkstemp(dir="test/test_rearr", suffix=".alg")

    def test_align(self):
        subprocess.run(
            f'''rearrangement < {self.post_file} 3< {self.ref_file} -s0 -6 -s1 4 -s2 2 -u -3 -v -9 -ru 0 -rv 0 -qu 0 -qv -5 > {self.alg_file}''',
            shell=True,
            executable="/bin/bash"
        )
        self.assertTrue(
            filecmp.cmp(
                self.alg_file,
                'test/test_rearr/output.alg',
                shallow=False
            )
        )
        subprocess.run(
            f'''gawk -f correct_micro_homology.awk -- {self.ref_file} {self.correct_file} < test/test_rearr/output.alg > {self.alg_file}.corrected''',
            shell=True,
            executable="/bin/bash"
        )
        self.assertTrue(
            filecmp.cmp(
                f'{self.alg_file}.corrected',
                'test/test_rearr/output.alg.corrected',
                shallow=False
            )
        )

    def tearDown(self):
        os.remove(self.ref_file)
        os.remove(self.correct_file)
        os.remove(self.post_file)
        os.remove(self.alg_file)
        os.remove(f'{self.alg_file}.corrected')

if __name__ == '__main__':
    unittest.main()
