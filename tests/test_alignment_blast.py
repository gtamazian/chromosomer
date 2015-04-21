#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import unittest
from chromosomer.alignment.blast import Blast
from chromosomer.alignment.exception import BlastAlignmentError

path = os.path.dirname(__file__)
os.chdir(path)


class TestBlastAlignment(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'blast', 'blast_tabular_alignments.txt'
        )
        self.__test_line = os.path.join(
            'data', 'blast', 'blast_tabular_line.txt'
        )
        self.__incorrect_lines = os.listdir(os.path.join(
            'data', 'blast', 'incorrect_input'
        ))
        # silence the logging messages
        logging.disable(logging.ERROR)


    def test_alignments(self):
        """
        Check if the parser reads a file in the BLAST tabular format
        in the correct way.
        """
        # test against the correct input file
        parser = Blast(self.__correct_file)
        for alignment in parser.alignments():
            self.assertEqual(len(alignment), 12)

        # test against the correct input line and check the entries
        parser = Blast(self.__correct_file)
        alignment = parser.alignments().next()
        self.assertIsInstance(alignment, Blast.Alignment)
        self.assertEqual(alignment.query, 'lcl|BA000007.2_gene_5374')
        self.assertEqual(alignment.subject,
                         'gi|556503834|ref|NC_000913.3|')
        self.assertEqual(alignment.identity, 98.14)
        self.assertEqual(alignment.length, 2856)
        self.assertEqual(alignment.mismatches, 53)
        self.assertEqual(alignment.gap_openings, 0)
        self.assertEqual(alignment.q_start, 1)
        self.assertEqual(alignment.q_end, 2856)
        self.assertEqual(alignment.s_start, 4483837)
        self.assertEqual(alignment.s_end, 4480982)
        self.assertEqual(alignment.e_value, 0.0)
        self.assertEqual(alignment.bit_score, 4981)

        # test again incorrect input lines, the exception must be
        # raised
        for incorrect_input in self.__incorrect_lines:
            parser = Blast(os.path.join(
                'data', 'blast', 'incorrect_input', incorrect_input))
            with self.assertRaises(BlastAlignmentError):
                alignment = parser.alignments().next()
                self.assertIsInstance(alignment, Blast.Alignment)

suite = unittest.TestLoader().loadTestsFromTestCase(TestBlastAlignment)
unittest.TextTestRunner(verbosity=2).run(suite)
