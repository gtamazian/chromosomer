#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import glob
import os
import tempfile
import unittest
from chromosomer.fasta import RandomSequence
from chromosomer.fasta import Writer
from chromosomer.wrapper.blast import MakeBlastDb

path = os.path.dirname(__file__)
os.chdir(path)


class TestWrapperMakeBlastDb(unittest.TestCase):
    def setUp(self):
        self.__fasta = tempfile.mkstemp()[1]
        self.__seq_length = 100
        self.__seq_number = 10

    def tearDown(self):
        temp_files = glob.glob('{}*'.format(self.__fasta))
        for i in temp_files:
            os.unlink(i)

    def test_launch(self):
        """
        The the makeblastdb launching routine.
        """
        with Writer(self.__fasta) as fasta_writer:
            for i in xrange(self.__seq_number):
                seq_generator = RandomSequence(self.__seq_length)
                fasta_writer.write('seq{}'.format(i+1),
                                   seq_generator.get())

        wrapper = MakeBlastDb(self.__fasta)
        wrapper.launch()

suite = unittest.TestLoader().loadTestsFromTestCase(
    TestWrapperMakeBlastDb)
unittest.TextTestRunner(verbosity=2).run(suite)
