#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import unittest
from chromosomer.alignment.lav import Lav

path = os.path.dirname(__file__)
os.chdir(path)


class TestLavAlignment(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'lav', 'lav_alignments.txt'
        )
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_alignments(self):
        """
        Check if the parser reads a file in the BLAST tabular format
        in the correct way.
        """
        # test against the correct input file
        parser = Lav(self.__correct_file)
        for alignment in parser.alignments():
            self.assertEqual(len(alignment), 6)

suite = unittest.TestLoader().loadTestsFromTestCase(TestLavAlignment)
unittest.TextTestRunner(verbosity=2).run(suite)
