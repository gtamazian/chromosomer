#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import pyfaidx
from chromosomer.exception import FragmentMapError
from collections import defaultdict
from collections import namedtuple
from operator import attrgetter

logging.basicConfig()
logger = logging.getLogger(__name__)


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
        with open(filename) as input_map_file:
            for line in input_map_file:
                line_parts = line.split('\t', 8)
                if len(line_parts) < 8:
                    logger.error('the incorrect number of columns')
                    raise FragmentMapError
                for i in self.numeric_values:
                    try:
                        line_parts[i] = int(line_parts[i])
                    except ValueError:
                        logger.error('the incorrect numeric value '
                                     '%s', line_parts[i])
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
