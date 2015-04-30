#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import unittest
from chromosomer.alignment.exception import LavAlignmentError
from chromosomer.alignment.lav import Lav

path = os.path.dirname(__file__)
os.chdir(path)


class TestLavAlignment(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'lav', 'lav_alignments.txt'
        )
        self.__incorrect_file_dir = os.path.join(
            'data', 'lav', 'incorrect_input'
        )
        self.__incorrect_files = os.listdir(self.__incorrect_file_dir)
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_alignments(self):
        """
        Check if the parser reads a file in the LAV format in the
        correct way.
        """
        # test against the correct input file
        parser = Lav(self.__correct_file)
        for alignment in parser.alignments():
            self.assertEqual(len(alignment), 6)
        # test againts incorrect input files
        for lav_file in self.__incorrect_files:
            parser = Lav(os.path.join(self.__incorrect_file_dir,
                                      lav_file))
            with self.assertRaises(LavAlignmentError):
                parser.alignments().next()

suite = unittest.TestLoader().loadTestsFromTestCase(TestLavAlignment)
unittest.TextTestRunner(verbosity=2).run(suite)
