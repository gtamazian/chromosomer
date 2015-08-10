#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import bioformats.bed
import bioformats.gff3
import csv
import logging
import vcf
from chromosomer.fragment import AlignmentToMap
from chromosomer.fragment import SeqLengths
from chromosomer.fragment import Map
from chromosomer.transfer import BedTransfer
from chromosomer.transfer import Gff3Transfer
from chromosomer.transfer import VcfTransfer
from bioformats.blast import BlastTab

from chromosomer.fragment import logger
logger.setLevel(logging.INFO)
logger.propagate = False
formatter = logging.Formatter('%(asctime)-15s - %(message)s',
                              '%Y-%m-%d %H:%M:%S')
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(formatter)
logger.addHandler(ch)

logging.basicConfig()
cli_logger = logging.getLogger(__name__)
cli_logger.propagate = False
cli_logger.addHandler(ch)
cli_logger.setLevel(logging.INFO)


def read_fragment_lengths(filename):
    """
    Given a name of a file with fragment lengths, read them to a
    dictionary.

    :param filename: a name of a file with fragment lengths
    :type filename: str
    :return: a dictionary which keys are fragment sequence names and
        values are their lengths
    :rtype: dict
    """
    result = dict()
    with open(filename) as length_file:
        length_reader = csv.reader(length_file, delimiter='\t')
        for fragment, length in length_reader:
            result[fragment] = int(length)
    return result


def chromosomer():
    """
    The main function that is run if Chromosomer was launched. It
    defines a command-line parser which processed arguments passed to
    the program.
    """
    parser = argparse.ArgumentParser(
        description='Reference-assisted chromosome assembly tool.')
    subparsers = parser.add_subparsers(dest='command')

    # Parser for the 'chromosomer assemble' part that produces a FASTA
    # file of assembled chromosomes from the specified fragment map.
    assemble_parser = subparsers.add_parser(
        'assemble',
        help='get sequences of assembled chromosomes',
        description='Get the FASTA file of assembled chromosomes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # required arguments for the 'assemble' routine
    assemble_parser.add_argument('map',
                                 help='a fragment map file')
    assemble_parser.add_argument('fragment_fasta',
                                 help='a FASTA file of fragment '
                                      'sequences to be assembled')
    assemble_parser.add_argument('output_fasta',
                                 help='the output FASTA file of the '
                                      'assembled chromosome sequences')

    # optinal arguments for the 'assemble' routine
    assemble_parser.add_argument('-s', '--save_soft_mask',
                                 action='store_true',
                                 help='keep soft masking from the '
                                      'original fragment sequences')

    # Parser for the 'chromosomer fragmentmap' part that
    # produces a map of fragment positions on reference
    # chromosomes from BLAST alignments of the fragments to the
    # chromosomes.
    fragmentmap_parser = subparsers.add_parser(
        'fragmentmap',
        description='Construct a fragment map from fragment '
                    'alignments to reference chromosomes.',
        help='construct a fragment map from alignments',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # required arguments for the 'fragmentmap' routine
    fragmentmap_parser.add_argument(
        'alignment_file',
        help='a BLAST tabular file of fragment alignments to '
             'reference chromosomes'
    )
    fragmentmap_parser.add_argument(
        'gap_size', type=int,
        help='a size of a gap inserted between mapped fragments'
    )
    fragmentmap_parser.add_argument(
        'fragment_lengths',
        help='a file containing lengths of fragment sequences; it can'
             'be obtained using \'chromosomer fastalength\''
    )
    fragmentmap_parser.add_argument(
        'output_map',
        help='an output fragment map file name'
    )

    # optional arguments for the 'fragmentmap' routine
    fragmentmap_parser.add_argument(
        '-r', '--ratio_threshold', type=float, default=1.2,
        help='the least ratio of two greatest fragment alignment '
             'scores to determine the fragment placed to a reference '
             'genome'
    )

    # Parser for the 'chromosomer transfer' part that transfers
    # genome feature annotation from fragments to their assembly
    transfer_parser = subparsers.add_parser(
        'transfer',
        description='Transfer annotated genomic features from '
                    'fragments to their assembly.',
        help='transfer annotated features from fragments to '
             'chromosomes',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # required arguments for the 'transfer' routine
    transfer_parser.add_argument('map',
                                 help='a fragment map file')
    transfer_parser.add_argument('annotation',
                                 help='a file of annotated genome '
                                      'features')
    transfer_parser.add_argument('output',
                                 help='an output file of the '
                                      'transfered annotation')

    # optional arguments for the 'transfer' routine
    transfer_parser.add_argument('-f', '--format', default='bed',
                                 help='the format of a file of '
                                      'annotated features (bed, '
                                      'gff3 or vcf)')

    # Parser for the 'chromosomer fastalength' part that calculates
    # lengths of sequences in the given FASTA file.
    fastalength_parser = subparsers.add_parser(
        'fastalength',
        description='Get lengths of sequences in the specified FASTA '
                    'file (required to build a fragment map)',
        help='get lengths of sequences from a FASTA file',
    )

    # required arguments for the 'fastalength' routine
    fastalength_parser.add_argument('fasta',
                                    help='a FASTA file which sequence '
                                         'lengths are to be obtained')
    fastalength_parser.add_argument('output',
                                    help='an output file of sequence '
                                         'lengths')

    args = parser.parse_args()

    if args.command == 'assemble':
        fragment_map = Map()
        fragment_map.read(args.map)
        fragment_map.assemble(args.fragment_fasta,
                              args.output_fasta,
                              args.save_soft_mask)
    elif args.command == 'fragmentmap':
        fragment_lengths = read_fragment_lengths(args.fragment_lengths)
        map_creator = AlignmentToMap(args.gap_size, fragment_lengths)
        alignments = BlastTab(args.alignment_file)
        fragment_map = map_creator.blast(alignments,
                                         args.ratio_threshold)
        fragment_map.write(args.output_map)
    elif args.command == 'transfer':
        total_count = transferred_count = 0
        if args.format == 'bed':
            transferrer = BedTransfer(args.map)
            with bioformats.bed.Writer(args.output) as output_file:
                for feature in bioformats.bed.Reader(
                        args.annotation).records():
                    total_count += 1
                    transferred_feature = transferrer.feature(feature)
                    if transferred_feature is not None:
                        transferred_count += 1
                        output_file.write(transferred_feature)
        elif args.format == 'gff3':
            transferrer = Gff3Transfer(args.map)
            with bioformats.gff3.Writer(args.output) as output_file:
                for feature in bioformats.gff3.Reader(
                        args.annotation).records():
                    total_count += 1
                    transferred_feature = transferrer.feature(feature)
                    if transferred_feature is not None:
                        transferred_count += 1
                        output_file.write(transferred_feature)
        elif args.format == 'vcf':
            transferrer = VcfTransfer(args.map)
            reader = vcf.Reader(open(args.annotation))
            writer = vcf.Writer(open(args.output, 'w'), reader)
            for variant in reader:
                total_count += 1
                transferred_feature = transferrer.feature(variant)
                if transferred_feature is not None:
                    transferred_count += 1
                    writer.write_record(transferred_feature)
            writer.close()

        cli_logger.info('%d features transferred', transferred_count)
        cli_logger.info('%d features skipped',
                        total_count - transferred_count)
    elif args.command == 'fastalength':
        seq_lengths = SeqLengths(args.fasta)
        with open(args.output, 'wt') as length_file:
            length_writer = csv.writer(length_file, delimiter='\t')
            for header, length in seq_lengths.lengths().iteritems():
                length_writer.writerow((header, length, ))
