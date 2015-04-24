#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import pyfaidx
import unittest
from chromosomer.fragmentmap import FastaWriter

path = os.path.dirname(__file__)
os.chdir(path)


class TestFastaWriter(unittest.TestCase):
    def setUp(self):
        self.__test_dir = os.path.join('data')

    def test_write(self):
        """
        Check if sequences are written to a FASTA file properly.
        """
        output_file = os.path.join(self.__test_dir, 'test.fa')

        # prepare the sequences
        fasta_patterns = ['AC', 'ACG', 'CCGT']
        sequences = {}
        for i, pattern in enumerate(fasta_patterns):
            sequences['chr{}'.format(i)] = pattern * 100

        with FastaWriter(output_file) as output_fasta:
            for header, sequence in sequences.iteritems():
                output_fasta.write(header, sequence)

        # check the written file
        reader = pyfaidx.Fasta(output_file)
        for header, sequence in sequences.iteritems():
            self.assertEqual(sequence, reader[header][:].seq)

        os.unlink(output_file)
        os.unlink(output_file + '.fai')

suite = unittest.TestLoader().loadTestsFromTestCase(TestFastaWriter)
unittest.TextTestRunner(verbosity=2).run(suite)
