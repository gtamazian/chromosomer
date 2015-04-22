#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
from collections import namedtuple
from exceptions import BedError

logging.basicConfig()
logger = logging.getLogger(__name__)

bed_columns = ('seq', 'start', 'end', 'name', 'score', 'strand',
               'thick_start', 'thick_end', 'color', 'block_num',
               'block_sizes', 'block_starts', 'extra')

bed_numeric_fields = (1, 2, 4, 6, 7, 9)

BedRecord = namedtuple('BedRecord', bed_columns)


class Reader(object):
    """
    This class implements a parser to read data from a file in the
    BED format.
    """
    def __init__(self, filename):
        """
        Given a name of a file, create a BED reader object to read
        data from it.

        :param filename: a name of a BED file
        :type filename: str
        """
        self.__filename = filename
        self.__lineno = 0
        self.__line = ''

    def records(self):
        """
        Iterate through records in the BED file the object was
        created from.

        :return: a record from the BED file the object was created from
        :rtype: BedRecord
        """
        with open(self.__filename) as bed_file:
            for self.__line in bed_file:
                self.__lineno += 1
                yield self.__parse_bed_line()

    def __parse_bed_line(self):
        """
        Parse the current line from the BED file.

        :return: a record from the BED file the object was created from
        :rtype: BedRecord
        """
        line_parts = self.__line.split()

        # convert numeric values: start, end, score, thick_start,
        # thick_end, block_num, blocks_sizes and block_starts; the
        # given BED line may contain the lesser number of columns,
        # so we adjust the tuple of numeric value positions
        line_numeric_pos = [x for x in bed_numeric_fields
                            if x < len(line_parts)]
        for i in line_numeric_pos:
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logger.error('line {}: the incorrect numeric value {'
                             '}'.format(self.__lineno, line_parts[i]))
                raise BedError

        # form the tuple to be returned as a result
        line_parts += ([None] * (len(bed_columns) - len(line_parts)))
        result = BedRecord(*line_parts)

        return result
