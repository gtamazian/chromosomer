#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import unittest
from chromosomer.track.bed import BedRecord
from chromosomer.track.bed import Reader

path = os.path.dirname(__file__)
os.chdir(path)


class TestBedReader(unittest.TestCase):
    def setUp(self):
        self.__correct_file = os.path.join(
            'data', 'bed', 'correct.bed'
        )
        # silence the logging messages
        logging.disable(logging.ERROR)

    def test_records(self):
        """
        Check if the parser reads a file in the BED format in the
        correct way.
        """
        # test against the correct input file
        parser = Reader(self.__correct_file)
        for record in parser.records():
            self.assertIsInstance(record, BedRecord)

suite = unittest.TestLoader().loadTestsFromTestCase(TestBedReader)
unittest.TextTestRunner(verbosity=2).run(suite)
