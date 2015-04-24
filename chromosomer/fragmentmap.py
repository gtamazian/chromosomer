#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import pyfaidx
import string
from chromosomer.exception import FragmentMapError
from collections import defaultdict
from collections import namedtuple
from operator import attrgetter

logging.basicConfig()
logger = logging.getLogger(__name__)


class FastaWriter(object):
    """
    The class implements routines to write sequences in the FASTA
    format.
    """

    def __init__(self, filename, width=72):
        """
        Create a FastaWriter object to write sequences in a FASTA file.

        :param filename: a name of a file to write sequences to
        :type filename: str
        """
        self.__filename = filename
        self.__width = width

    def __enter__(self):
        self.__output = open(self.__filename, 'w')
        return self

    def write(self, header, sequence):
        """
        Write a sequence to the FASTA file.

        :param header: a sequence header
        :param sequence: a sequence
        :type header: str
        :type sequence: str
        """
        self.__output.write('>{}\n'.format(header))
        seq_lines = []
        for i in xrange(0, len(sequence), self.__width):
            seq_lines.append(sequence[i:i + self.__width])
        self.__output.write('\n'.join(seq_lines))
        self.__output.write('\n')

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__output.close()


class FragmentMap(object):
    """
    The class implements routines related to creation, reading and
    writing a fragment map that describes how genome fragments are
    situated on reference genome chromosomes.
    """

    numeric_values = (1, 2, 3, 6, 7)
    
    record_names = ('fr_name', 'fr_length', 'fr_start', 'fr_end',
                    'fr_strand', 'ref_chr', 'ref_start', 'ref_end')

    Record = namedtuple('Record', record_names)

    def __init__(self):
        """
        Initializes a FragmentMap object.
        """
        self.__fragments = defaultdict(list)
        self.__block_adding = False

    def add_record(self, new_record):
        """
        Given a new fragment record, add it to the fragment map.

        :param new_record: a record to be added to the map
        :type new_record: FragmentMap.Record
        """
        self.__fragments[new_record.ref_chr].append(new_record)

    def read(self, filename):
        """
        Read a fragment map from the specified file. The file records
        are added to the records in the map.

        :param filename: a name of a file to read a fragment map from
        :type: str
        """
        lineno = 0
        with open(filename) as input_map_file:
            for line in input_map_file:
                lineno += 1
                line_parts = line.split('\t', 8)
                if len(line_parts) < 8:
                    logger.error('line %d: the incorrect number of '
                                 'columns', lineno)
                    raise FragmentMapError
                for i in self.numeric_values:
                    try:
                        line_parts[i] = int(line_parts[i])
                    except ValueError:
                        logger.error('line %d: the incorrect numeric '
                                     'value %s', lineno, line_parts[i])
                        raise FragmentMapError
                new_record = FragmentMap.Record(*line_parts)
                self.add_record(new_record)

    def chromosomes(self):
        """
        Return an iterator to the chromosomes the fragment map
        describes.

        :return: an iterator to iterate through the fragment map
            chromosomes
        """
        sorted_chromosomes = sorted(self.__fragments.keys())
        for i in sorted_chromosomes:
            yield i

    def fragments(self, chromosome):
        """
        Return an iterator to fragments of the specified chromosome
        describes by the map. If the chromosome is absent, None is
        returned.

        :param chromosome: a chromosome which fragments are to be
            iterated
        :type chromosome: str
        :return: an iterator to iterate through the chromosome's
            fragments
        """
        if chromosome not in self.__fragments:
            logging.error('%s missing in the fragment map', chromosome)
            raise FragmentMapError

        sorted_fragments = sorted(self.__fragments[chromosome],
                                  key=attrgetter('ref_start'))

        for i in sorted_fragments:
            yield i

    def write(self, filename):
        """
        Write the fragment map to the specified file.

        :param filename: a name of a file to write the fragment map to
        :type filename: str
        """
        template = '\t'.join(['{}'] * len(self.record_names)) + '\n'
        with open(filename, 'w') as output_map_file:
            for chromosome in self.chromosomes():
                for fragment in self.fragments(chromosome):
                    new_line = template.format(*fragment)
                    output_map_file.write(new_line)

    def assebmle(self, fragment_filename, output_filename):
        """
        Assemble chromosome sequences from fragments.

        :param fragment_filename: a name of a FASTA file of fragment
            sequences
        :param output_filename: a name of the output FASTA file of
            the assembled chromosomes
        """
        fragment_fasta = pyfaidx.Fasta(fragment_filename)
        complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
        with FastaWriter(output_filename) as chromosome_writer:
            for chromosome in self.chromosomes():
                seq = []
                for record in self.fragments(chromosome):
                    if record.fr_name == 'GAP':
                        record_seq = 'N' * (record.fr_end -
                                            record.fr_start)
                    else:
                        if record.fr_name not in fragment_fasta:
                            logger.error('the fragment %s sequence '
                                         'missing', record.fr_name)
                            raise FragmentMapError
                        record_seq = fragment_fasta[record.fr_name][
                            record.fr_start:record.fr_end].seq
                        # convert the sequence to non-unicode
                        record_seq = str(record_seq)
                        # if the fragment orientation is reverse, then
                        # the reverse complement of the fragment
                        # sequence is written
                        if record.fr_strand == '-':
                            record_seq = record_seq.translate(
                                complement)
                    seq.append(record_seq)
                chromosome_writer.write(chromosome, ''.join(seq))
