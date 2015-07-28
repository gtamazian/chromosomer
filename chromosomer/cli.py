#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
from chromosomer.fragment import AlignmentToMap
from chromosomer.fragment import SeqLengths
from chromosomer.fragment import Map
from bioformats.blast import BlastTab


def chromosomer():
    """
    The main function that is run if Chromosomer was launched. It
    defines a command-line parser which processed arguments passed to
    the program.
    """
    parser = argparse.ArgumentParser(description='Reference-assisted '
                                                 'chromosome assembly '
                                                 'tool.')
    subparsers = parser.add_subparsers(dest='command')

    # Parser for the 'chromosomer assemble' part that produces a FASTA
    # file of assembled chromosomes from the specified fragment map.
    assemble_parser = subparsers.add_parser(
        'assemble',
        description='Get the FASTA file of assembled chromosomes.'
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

    # required arguments for the 'chromosomer fragmentmap' part that
    # produces a map of fragment positions on reference
    # chromosomes from BLAST alignments of the fragments to the
    # chromosomes.
    fragmentmap_parser = subparsers.add_parser(
        'fragmentmap',
        description='Construct a fragment map from fragment '
                    'alignments to reference chromosomes'
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
        help='a file containing lengths of fragment sequences'
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

    args = parser.parse_args()

    if args.command == 'assemble':
        fragment_map = Map()
        fragment_map.read(args.map)
        fragment_map.assemble(args.fragment_fasta,
                              args.output_fasta)
    elif args.command == 'fragmentmap':
        fragment_lengths = SeqLengths(args.fragment_lengths)
        map_creator = AlignmentToMap(args.gap_size,
                                     fragment_lengths.lengths())
        alignments = BlastTab(args.alignment_file)
        fragment_map = map_creator.blast(alignments,
                                         args.ratio_threshold)
        fragment_map.write(args.output_map)
