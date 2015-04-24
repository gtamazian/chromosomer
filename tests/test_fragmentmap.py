#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import unittest
from chromosomer.fragmentmap import FragmentMap
from itertools import izip

path = os.path.dirname(__file__)
os.chdir(path)


class TestFragmentMap(unittest.TestCase):
    def setUp(self):
        self.__test_line = os.path.join(
            'data', 'fragment_map', 'fragment_map_line.txt'
        )
        self.__output_file = os.path.join(
            'data', 'fragment_map', 'fragment_map_output.txt'
        )
        # silence the logging message
        logging.disable(logging.ERROR)

    def test_add_record(self):
        """
        Check if fragment records are added correctly.
        """
        fragment_map = FragmentMap()
        new_record = FragmentMap.Record(
            fr_name='fragment1',
            fr_length=180,
            fr_start=0,
            fr_end=180,
            fr_strand='+',
            ref_chr='chr1',
            ref_start=5000,
            ref_end=5180
        )
        fragment_map.add_record(new_record)

    def test_read(self):
        """
        Test the FragmentMap reading routine.
        """
        fragment_map = FragmentMap()
        fragment_map.read(self.__test_line)
        fragment = fragment_map.fragments('chr1').next()
        self.assertEqual(fragment.fr_name, 'fragment1')
        self.assertEqual(fragment.fr_length, 180)
        self.assertEqual(fragment.fr_start, 0)
        self.assertEqual(fragment.fr_end, 180)
        self.assertEqual(fragment.fr_strand, '+')
        self.assertEqual(fragment.ref_chr, 'chr1')
        self.assertEqual(fragment.ref_start, 5000)
        self.assertEqual(fragment.ref_end, 5180)

    def test_chromosomes(self):
        """
        Test the FragmentMap chromosomes iterator.
        """
        fragment_map = FragmentMap()
        fragment_map.read(self.__test_line)
        chromosomes = list(fragment_map.chromosomes())
        self.assertEqual(chromosomes, ['chr1'])

    def test_fragments(self):
        """
        Test the FragmentMap fragments iterator.
        """
        fragment_map = FragmentMap()
        fragment_map.read(self.__test_line)
        fragments = list(fragment_map.fragments('chr1'))
        self.assertEqual(len(fragments), 1)
        self.assertIsInstance(fragments[0], FragmentMap.Record)

    def test_write(self):
        """
        Test the FragmentMap writing routine.
        """
        fragment_map = FragmentMap()
        fragment_map.read(self.__test_line)
        fragment_map.write(self.__output_file)

        with open(self.__output_file) as output_file:
            with open(self.__test_line) as original_file:
                for x, y in izip(original_file, output_file):
                    self.assertEqual(x, y)

        os.unlink(self.__output_file)

suite = unittest.TestLoader().loadTestsFromTestCase(TestFragmentMap)
unittest.TextTestRunner(verbosity=2).run(suite)
