#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import os
import glob
import logging
import pyfaidx
import string
import tempfile
import unittest
from chromosomer.alignment.blast import Blast
from chromosomer.fasta import RandomSequence
from chromosomer.fasta import Writer
from chromosomer.fragment import AlignmentToMap
from chromosomer.fragment import Length
from chromosomer.fragment import Map
from chromosomer.fragment import MapError
from chromosomer.fragment import Simulator
from chromosomer.wrapper.blast import BlastN
from chromosomer.wrapper.blast import MakeBlastDb
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
        fragment_map = Map()
        new_record = Map.Record(
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
        Test the Map reading routine.
        """
        fragment_map = Map()
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
            with self.assertRaises(MapError):
                fragment_map.read(os.path.join(
                    self.__incorrect_file_dir, i))

    def test_chromosomes(self):
        """
        Test the Map chromosomes iterator.
        """
        fragment_map = Map()
        fragment_map.read(self.__test_line)
        chromosomes = list(fragment_map.chromosomes())
        self.assertEqual(chromosomes, ['chr1'])

    def test_fragments(self):
        """
        Test the Map fragments iterator.
        """
        fragment_map = Map()
        fragment_map.read(self.__test_line)
        fragments = list(fragment_map.fragments('chr1'))
        self.assertEqual(len(fragments), 1)
        self.assertIsInstance(fragments[0], Map.Record)

        # check if the missing chromosome is processed correctly
        with self.assertRaises(MapError):
            list(fragment_map.fragments('chrN'))

    def test_write(self):
        """
        Test the Map writing routine.
        """
        fragment_map = Map()
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
                    chromosomes[i].append(fr_seq[::-1].translate(
                        complement))
                else:
                    chromosomes[i].append(fr_seq)
                chromosomes[i].append('N' * gap_size)
            chromosomes[i] = ''.join(chromosomes[i])
        # contruct a fragment __map
        fragment_map = Map()
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
                fragment_map.add_record(Map.Record(
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
                fragment_map.add_record(Map.Record(
                    fr_name, fr_length, fr_start, fr_end, fr_strand,
                    ref_chr, ref_start, ref_end
                ))
                current_start += fr_length

        output_chromosomes = os.path.join(self.__output_dir,
                                          'temp_chromosomes.txt')
        output_fragments = os.path.join(self.__output_dir,
                                        'temp_fragments.txt')

        # write the fragment sequences to a FASTA file
        with Writer(output_fragments) as writer:
            for i, j in fragments.iteritems():
                writer.write(i, j)

        fragment_map.assemble(output_fragments, output_chromosomes)

        # read fragments from the written FASTA file and compare them
        # to the original ones
        assembled_chromosomes = pyfaidx.Fasta(output_chromosomes)
        for i, seq in chromosomes.iteritems():
            self.assertEqual(seq, assembled_chromosomes[i][:].seq)

        # try to use the fragment absent in the FASTA file of
        # fragment sequences
        fragment_map.add_record(Map.Record(
            fr_name='missing_fragment',
            fr_length=0,
            fr_start=0,
            fr_end=0,
            fr_strand='+',
            ref_chr='chr3',
            ref_start=0,
            ref_end=0
        ))
        with self.assertRaises(MapError):
            fragment_map.assemble(output_fragments, output_chromosomes)

        os.unlink(output_chromosomes)
        os.unlink(output_chromosomes + '.fai')
        os.unlink(output_fragments)
        os.unlink(output_fragments + '.fai')


class TestFragmentLength(unittest.TestCase):
    def setUp(self):
        self.__fragment_number = 10
        self.__fragment_length = 10
        self.__fasta_temp = tempfile.mkstemp()[1]

        # create a FASTA file of random sequences
        with Writer(self.__fasta_temp) as fasta_writer:
            seq_generator = RandomSequence(self.__fragment_length)
            for i in xrange(self.__fragment_length):
                fasta_writer.write('seq{}'.format(i+1),
                                   seq_generator.get())

    def test_lengths(self):
        """
        Test the lengths method.
        """
        x = Length(self.__fasta_temp)
        lengths = x.lengths()
        for i in lengths.itervalues():
            self.assertEqual(i, self.__fragment_length)

    def tearDown(self):
        os.unlink(self.__fasta_temp)


class TestFragmentSimulator(unittest.TestCase):
    def setUp(self):
        self.__fragment_number = 10
        self.__chromosome_number = 2
        self.__fragment_length = 10
        self.__gap_size = 5

        self.__simulator = Simulator(self.__fragment_length,
                                     self.__fragment_number,
                                     self.__chromosome_number,
                                     self.__gap_size)

    def test_write(self):
        """
        The the writing method of the fragment simulator.
        """
        self.__fragments = tempfile.mkstemp()[1]
        self.__chromosomes = tempfile.mkstemp()[1]
        self.__map = tempfile.mkstemp()[1]

        self.__simulator.write(self.__map, self.__fragments,
                               self.__chromosomes)

        # check if the correct number of fragment and chromosome
        # sequences was written
        fragment_fasta = pyfaidx.Fasta(self.__fragments)
        self.assertEqual(len(fragment_fasta.keys()),
                         self.__fragment_number)
        chromosome_fasta = pyfaidx.Fasta(self.__chromosomes)
        self.assertEqual(len(chromosome_fasta.keys()),
                         self.__chromosome_number)

        # check if a correct fragment map was written
        test_map = Map()
        test_map.read(self.__map)

        os.unlink(self.__fragments)
        os.unlink(self.__fragments + '.fai')
        os.unlink(self.__chromosomes)
        os.unlink(self.__chromosomes + '.fai')
        os.unlink(self.__map)


class TestFragmentAlignmentToMap(unittest.TestCase):
    def setUp(self):
        # simulate fragments and chromosomes
        self.__fragment_number = 10
        self.__chromosome_number = 2
        self.__fragment_length = 100
        self.__gap_size = 5

        self.__simulator = Simulator(self.__fragment_length,
                                     self.__fragment_number,
                                     self.__chromosome_number,
                                     self.__gap_size)

        # create the corresponding files
        self.__map_file = tempfile.mkstemp()[1]
        self.__fragment_file = tempfile.mkstemp()[1]
        self.__chromosome_file = tempfile.mkstemp()[1]

        self.__simulator.write(self.__map_file, self.__fragment_file,
                               self.__chromosome_file)

        # create the chromosome sequence database and align the
        # fragments to it
        makeblastdb_wrapper = MakeBlastDb(self.__chromosome_file)
        makeblastdb_wrapper.launch()

        self.__alignment_file = tempfile.mkstemp()[1]
        blastn_wrapper = BlastN(self.__fragment_file,
                                self.__chromosome_file,
                                self.__alignment_file)
        blastn_wrapper.set('-outfmt', 6)
        blastn_wrapper.set('-dust', 'no')
        blastn_wrapper.launch()

    def test_blast(self):
        """
        Test the blast method which utilizes BLASTN alignments to
        construct a fragment map.
        """
        fragment_lengths = Length(self.__fragment_file)
        map_creator = AlignmentToMap(self.__gap_size,
                                     fragment_lengths.lengths())
        blast_alignments = Blast(self.__alignment_file)
        new_map = map_creator.blast(blast_alignments, 1.2)
        orig_map = Map()
        orig_map.read(self.__map_file)

        # compare the obtained fragment map with the original one
        for chromosome in orig_map.chromosomes():
            for orig, new in izip(orig_map.fragments(chromosome),
                                  new_map.fragments(chromosome)):
                self.assertEqual(orig, new)

    def tearDown(self):
        os.unlink(self.__map_file)
        os.unlink(self.__fragment_file)
        os.unlink(self.__chromosome_file)
        os.unlink(self.__alignment_file)
        for i in glob.glob('{}*'.format(self.__chromosome_file)):
            os.unlink(i)


suite = unittest.TestLoader().loadTestsFromTestCase(
    TestFragmentAlignmentToMap)
unittest.TextTestRunner(verbosity=2).run(suite)

suite = unittest.TestLoader().loadTestsFromTestCase(TestFragmentLength)
unittest.TextTestRunner(verbosity=2).run(suite)

suite = unittest.TestLoader().loadTestsFromTestCase(TestFragmentMap)
unittest.TextTestRunner(verbosity=2).run(suite)

suite = unittest.TestLoader().loadTestsFromTestCase(
    TestFragmentSimulator)
unittest.TextTestRunner(verbosity=2).run(suite)
