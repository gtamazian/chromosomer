#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import logging
import pyfaidx
import string
import unittest
from chromosomer.fragmentmap import FastaWriter
from chromosomer.fragmentmap import FragmentMap
from chromosomer.fragmentmap import FragmentMapError
from itertools import izip

path = os.path.dirname(__file__)
os.chdir(path)


class TestFragmentMap(unittest.TestCase):
    def setUp(self):
        self.__test_line = os.path.join(
            'data', 'fragment_map', 'fragment_map_line.txt'
        )
        self.__incorrect_file_dir = os.path.join(
            'data', 'fragment_map', 'incorrect_input'
        )
        self.__output_dir = os.path.join(
            'data', 'fragment_map'
        )
        self.__incorrect_files = os.listdir(self.__incorrect_file_dir)
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

        # check for incorrect input files
        for i in self.__incorrect_files:
            with self.assertRaises(FragmentMapError):
                fragment_map.read(os.path.join(
                    self.__incorrect_file_dir, i))

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

        # check if the missing chromosome is processed correctly
        with self.assertRaises(FragmentMapError):
            list(fragment_map.fragments('chrN'))

    def test_write(self):
        """
        Test the FragmentMap writing routine.
        """
        fragment_map = FragmentMap()
        fragment_map.read(self.__test_line)

        output_filename = os.path.join('data', 'fragment_map',
                                       'fragment_map_output.txt')
        fragment_map.write(output_filename)

        with open(output_filename) as output_file:
            with open(self.__test_line) as original_file:
                for x, y in izip(original_file, output_file):
                    self.assertEqual(x, y)

        os.unlink(output_filename)

    def test_assemble(self):
        """
        Test the assemble routine.
        """
        # first, we form fragment and chromosome sequences
        fragments = {}
        fragment_pattern = ['AC', 'AG', 'CT', 'CG', 'AT']
        for i, pattern in enumerate(fragment_pattern):
            fragments['fragment{}'.format(i+1)] = pattern * 5
        # a negative number indicated reverse orientation of a fragment
        chromosome_content = {'chr1': [1, -2, 3], 'chr2': [-4, 5]}
        # get chromosome sequences
        chromosomes = {}
        complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
        gap_size = 10
        for i, chromosome_fragments in chromosome_content.iteritems():
            chromosomes[i] = []
            for j in chromosome_fragments:
                fr_seq = fragments['fragment{}'.format(abs(j))]
                if j < 0:
                    chromosomes[i].append(fr_seq.translate(complement))
                else:
                    chromosomes[i].append(fr_seq)
                chromosomes[i].append('N' * gap_size)
            chromosomes[i] = ''.join(chromosomes[i])
        # contruct a fragment map
        fragment_map = FragmentMap()
        for i, chromosome_fragments in chromosome_content.iteritems():
            current_start = 0
            for j in chromosome_fragments:
                fr_name = 'fragment{}'.format(abs(j))
                fr_length = 10
                fr_start = 0
                fr_end = fr_length
                fr_strand = '+' if j > 0 else '-'
                ref_chr = i
                ref_start = current_start
                ref_end = current_start + fr_length
                fragment_map.add_record(FragmentMap.Record(
                    fr_name, fr_length, fr_start, fr_end, fr_strand,
                    ref_chr, ref_start, ref_end
                ))
                current_start += fr_length
                # add the gap
                fr_name = 'GAP'
                fr_length = gap_size
                fr_start = 0
                fr_end = gap_size
                fr_strand = '+'
                ref_chr = i
                ref_start = current_start
                ref_end = current_start + fr_end
                fragment_map.add_record(FragmentMap.Record(
                    fr_name, fr_length, fr_start, fr_end, fr_strand,
                    ref_chr, ref_start, ref_end
                ))
                current_start += fr_length

        output_chromosomes = os.path.join(self.__output_dir,
                                          'temp_chromosomes.txt')
        output_fragments = os.path.join(self.__output_dir,
                                        'temp_fragments.txt')

        # write the fragment sequences to a FASTA file
        with FastaWriter(output_fragments) as writer:
            for i, j in fragments.iteritems():
                writer.write(i, j)

        fragment_map.assebmle(output_fragments, output_chromosomes)

        # read fragments from the written FASTA file and compare them
        # to the original ones
        assembled_chromosomes = pyfaidx.Fasta(output_chromosomes)
        for i, seq in chromosomes.iteritems():
            self.assertEqual(seq, assembled_chromosomes[i][:].seq)

        # try to use the fragment absent in the FASTA file of
        # fragment sequences
        fragment_map.add_record(FragmentMap.Record(
            fr_name='missing_fragment',
            fr_length=0,
            fr_start=0,
            fr_end=0,
            fr_strand='+',
            ref_chr='chr3',
            ref_start=0,
            ref_end=0
        ))
        with self.assertRaises(FragmentMapError):
            fragment_map.assebmle(output_fragments, output_chromosomes)

        os.unlink(output_chromosomes)
        os.unlink(output_chromosomes + '.fai')
        os.unlink(output_fragments)
        os.unlink(output_fragments + '.fai')


suite = unittest.TestLoader().loadTestsFromTestCase(TestFragmentMap)
unittest.TextTestRunner(verbosity=2).run(suite)
