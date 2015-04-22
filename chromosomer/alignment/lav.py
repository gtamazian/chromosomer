#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from collections import namedtuple
from exception import LavAlignmentError
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)


class Lav(object):
    """
    This class implements a parse to read alignment files in the LAV
    format.
    """
    stanza_prefixes = ('s', 'h', 'a', 'x', 'm', 'census')

    Alignment = namedtuple('Alignment', stanza_prefixes)

    SStanzaItem = namedtuple('SStanzaItem', (
        'filename', 'start', 'stop', 'rev_comp_flag', 'seq_number'
    ))

    HStanzaItem = namedtuple('HStanzaItem', (
        'seq_name', 'rev_comp_flag'
    ))

    AStanzaItem = namedtuple('AStanzaItem', (
        'score', 'first_start', 'first_end', 'second_start',
        'second_end', 'segments'
    ))

    GapFreeSegment = namedtuple('GapFreeSegment', (
        'first_start', 'second_start', 'first_end', 'second_end',
        'identity'
    ))

    MStanza = namedtuple('MStanza', ('regions', 'base_count'))

    MStanzaItem = namedtuple('MStanzaItem', ('start', 'end'))

    GapFreeAlignment = namedtuple('GapFreeAlignment', (
        'first_seq', 'second_seq', 'first_start', 'first_end',
        'second_start', 'second_end'))

    def __init__(self, filename):
        """
        Given a name of a file, create a LAV alignment format parser
        from it.

        :param filename: a name of a file in the LAV format
        :type filename: str
        """
        self.__filename = filename
        self.__handler = None
        self.__lineno = 0
        self.__line = None

    def alignments(self):
        """
        Iterate through alignments in the file the object was created
        from.

        :return: an alignment from the file the object was created from
        :rtype: Lav.Alignment
        """
        with open(self.__filename) as self.__handler:
            self.__line = self.__handler.readline().rstrip()
            self.__lineno += 1
            # first, we read the d-stanza
            self.__parse_section_header()
            self.__parse_d_stanza()
            # read alignments
            self.__line = self.__handler.readline().rstrip()
            self.__lineno += 1
            while self.__parse_section_header():
                self.__parse_s_stanza()
                h_stanza = self.__parse_h_stanza()
                first_seq = h_stanza[0].seq_name.split()[0]
                second_seq = h_stanza[1].seq_name.split()[0]
                self.__line = self.__handler.readline().rstrip()
                self.__lineno += 1
                # read a stanza of the corresponding type according
                # to the header that we read
                while not self.__line.startswith('#'):
                    if self.__line.startswith('a'):
                        a_stanza = self.__parse_a_stanza()
                        for segment in a_stanza.segments:
                            assert isinstance(segment,
                                              Lav.GapFreeSegment)
                            gap_free_alignment = Lav.GapFreeAlignment(
                                first_seq=first_seq,
                                second_seq=second_seq,
                                first_start=segment.first_start,
                                first_end=segment.first_end,
                                second_start=segment.second_start,
                                second_end=segment.second_end
                            )
                            yield gap_free_alignment
                    elif self.__line.startswith('x'):
                        self.__parse_x_stanza()
                    elif self.__line.startswith('m'):
                        self.__parse_m_stanza()
                    elif self.__line.startswith('Census'):
                        self.__parse_census_stanza()
                    self.__line = self.__handler.readline().rstrip()
                    self.__lineno += 1

    def __parse_section_header(self):
        """
        Parse a section header, return False if the header indicates
        the end of the file.

        :return: a boolean value indicating if the current header
            corresponds to a new section
        :rtype: bool
        """
        if self.__line == '#:lav':
            return True
        elif self.__line == '#:eof':
            return False
        else:
            logger.error('self.__line {}: an incorrect section '
                         'header'.format(self.__lineno))
            raise LavAlignmentError

    def __parse_d_stanza(self):
        """
        Parse the d-stanza from a LAV file.

        :return: a string representing the d-stanza contents
        :rtype: str
        """
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        if self.__line != 'd {':
            logger.error('line {}: an incorrect d-stanza '
                         'header'.format(self.__lineno))
            raise LavAlignmentError

        # start reading the d-stanza contents
        stanza_contents = []
        self.__line = self.__handler.readline()
        self.__lineno += 1
        while self.__line.rstrip() != '}':
            stanza_contents.append(self.__line)
            self.__line = self.__handler.readline()
            self.__lineno += 1

        return ''.join(stanza_contents)

    def __parse_s_stanza(self):
        """
        Parse the s-stanza from a LAV file.

        :return: a list of SStanzaItem tuples representing the
            s-stanza contents
        :rtype: tuple
        """
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1

        # check the first line of the stanza
        if self.__line != 's {':
            logger.error('line {}: an incorrect s-stanza '
                         'header'.format(self.__lineno))
            raise LavAlignmentError

        first_line = self.__lineno

        # start reading the s-stanza
        sequence_files = []
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        while self.__line != '}':
            # there must be only two sequence descriptions in the
            # stanza
            if len(sequence_files) > 2:
                logger.error('line {}: too many sequence files in'
                             'the s stanza'.format(self.__lineno))
                raise LavAlignmentError
            # parse each sequence line; 5 is the greatest number of
            # parts in the sequence description
            line_parts = self.__line.split(None, 5)
            # there must be either 3 or 5 parts in the line
            if len(line_parts) not in {3, 5}:
                logger.error('line {}: the incorrect number of '
                             'sequence file description parts ({'
                             '})'.format(self.__lineno,
                                         len(line_parts)))
                raise LavAlignmentError
            # remove quotes from the file name
            if len(line_parts) > 2:
                line_parts[0] = line_parts[0][1:-1]
            else:
                logger.error('line {}: an empty file name'.format(
                    self.__lineno))
                raise LavAlignmentError
            # convert numeric values
            try:
                line_parts[1] = int(line_parts[1])
                line_parts[2] = int(line_parts[2])
                if len(line_parts) > 3:
                    line_parts[3] = int(line_parts[3])
                    line_parts[4] = int(line_parts[4])
            except ValueError:
                logger.error('line {}: incorrect numeric'
                             'values'.format(self.__lineno))
                raise LavAlignmentError
            # form a SStanzaField tuple and write it to the sequence
            # list
            if len(line_parts) < 5:
                line_parts = line_parts + [None] * 2
            sequence_files.append(Lav.SStanzaItem(*line_parts))
            # read the next line
            self.__line = self.__handler.readline().rstrip()
            self.__lineno += 1

        # check that there are exactly two sequence files
        if len(sequence_files) < 2:
            logger.error('lines {}--{}: a pair of sequence lines '
                         'required'.format(first_line, self.__lineno))
            raise LavAlignmentError

        return sequence_files

    def __parse_h_stanza(self):
        """
        Parse the h-stanza from a LAV file.

        :return: a named tuple of the HStanzaItem type representing
            the h-stanza contents
        :rtype: Lav.HStanzaItem
        """
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1

        # check the first line of the stanza
        if self.__line != 'h {':
            logger.error('line {}: an incorrect h-stanza '
                         'header'.format(self.__lineno))
            raise LavAlignmentError

        first_line = self.__lineno

        # start reading the h stanza
        header_lines = []
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        while self.__line != '}':
            # there must be only two sequence descriptions in the
            # stanza
            if len(header_lines) > 2:
                logger.error('line {}: too many header lines in the '
                             'h stanza'.format(self.__lineno))
                raise LavAlignmentError
            # parse the header line having been read
            line_parts = self.__line.split()
            # remove quote symbols
            if len(line_parts[0]) > 2:
                line_parts[0] = line_parts[0][1:-1]
            else:
                logger.error('line {}: an empty sequence '
                             'header'.format(self.__lineno))
                raise LavAlignmentError
            # each header line must start with '>'
            if not line_parts[0].startswith('>'):
                logger.error('line {}: the header line must start '
                             'with >'.format(self.__lineno))
                raise LavAlignmentError
            # check if a sequence name is present in the line
            if len(line_parts[0]) < 2:
                logger.error('line {}: the header line contains no '
                             'sequence name'.format(self.__lineno))
                raise LavAlignmentError
            else:
                # remove '>' from the header
                line_parts[0] = line_parts[0][1:]
            # check if there is the reverse complement flag in the line
            if len(line_parts) > 2:
                reverse_flag = line_parts[-2] == '(reverse' and \
                    line_parts[-1] == 'complement)'
            else:
                reverse_flag = False
            # cut line parts with the '>' sign and the reverse
            # complement flag if specified
            if reverse_flag:
                line_parts = line_parts[:-2]

            header_lines.append(Lav.HStanzaItem(
                seq_name=' '.join(line_parts),
                rev_comp_flag=reverse_flag))

            # read the next line
            self.__line = self.__handler.readline().rstrip()
            self.__lineno += 1

        # check that there are exactly two header lines
        if len(header_lines) < 2:
            logger.error('lines {}--{}: a pair of header lines '
                         'required'.format(first_line + 1,
                                           self.__lineno))
            raise LavAlignmentError

        return header_lines

    def __parse_a_stanza(self):
        """
        Parse the a-stanza from a LAV file.

        :return: a named tuple of the AStanzaItem type representing
            the a-stanza contents
        :rtype: Lav.AStanzaItem
        """
        # check the first line of the a stanza
        if self.__line != 'a {':
            logger.error('line {}: an incorrect a-stanza '
                         'header'.format(self.__lineno))
            raise LavAlignmentError

        # start reading the a stanza: read three lines describing the
        # alignment block

        # score
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        line_parts = self.__line.split(None, 2)
        if len(line_parts) == 2 and line_parts[0] == 's':
            try:
                alignment_block_score = int(line_parts[1])
            except ValueError:
                logger.error('line {}: the incorrect score value {'
                             '}'.format(self.__lineno, line_parts[1]))
                raise LavAlignmentError
        else:
            logger.error('line {}: the incorrect score record in an a-'
                         'stanza'.format(self.__lineno))
            raise LavAlignmentError

        # start positions of aligned sequences
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        line_parts = self.__line.split(None, 3)
        if len(line_parts) == 3 and line_parts[0] == 'b':
            try:
                alignment_block_first_start = int(line_parts[1])
                alignment_block_second_start = int(line_parts[2])
            except ValueError:
                logger.error('line {}: the incorrect position value {'
                             '} or {}'.format(self.__lineno,
                                              *line_parts))
                raise LavAlignmentError
        else:
            logger.error('line {}: an incorrect aligned block record '
                         'in an a-stanza'.format(self.__lineno))
            raise LavAlignmentError

        # end positions of aligned sequences
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        line_parts = self.__line.split(None, 3)
        if len(line_parts) == 3 and line_parts[0] == 'e':
            try:
                alignment_block_first_end = int(line_parts[1])
                alignment_block_second_end = int(line_parts[2])
            except ValueError:
                logger.error('line {}: an incorrect position value {'
                             '} or {}'.format(self.__lineno,
                                              *line_parts))
                raise LavAlignmentError
        else:
            logger.error('line {}: an incorrect aligned block record '
                         'in an a-stanza'.format(self.__lineno))
            raise LavAlignmentError

        alignment_block = Lav.AStanzaItem(
            score=alignment_block_score,
            first_start=alignment_block_first_start,
            first_end=alignment_block_first_end,
            second_start=alignment_block_second_start,
            second_end=alignment_block_second_end,
            segments=[])
        # read gap-free segments of the alignment block
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        while self.__line != '}':
            line_parts = self.__line.split(None, 6)
            if len(line_parts) == 6 and line_parts[0] == 'l':
                for i in xrange(1, 6):
                    try:
                        line_parts[i] = int(line_parts[i])
                    except ValueError:
                        logger.error(
                            'line {}: the incorrect numeric value {'
                            '}'.format(self.__lineno, line_parts[i]))
                        raise LavAlignmentError
                alignment_block.segments.append(
                    Lav.GapFreeSegment(*line_parts[1:]))
            else:
                logger.error('line {}: the incorrect segment record '
                             'in an a-stanza'.format(self.__lineno))
                raise LavAlignmentError
            # read the next line
            self.__line = self.__handler.readline().rstrip()
            self.__lineno += 1

        return alignment_block

    def __parse_x_stanza(self):
        """
        Parse the x-stanza from a LAV file.

        :return: the number of newly masked bases from the x-stanza
        :rtype: int
        """
        if self.__line != 'x {':
            logger.error('line {0}: an incorrect x-stanza header')
            raise LavAlignmentError

        # read a sigle line containing the number of newly masked bases
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        line_parts = self.__line.split(None, 2)
        newly_masked_bases = None
        if len(line_parts) == 2 and line_parts[0] == 'n':
            try:
                newly_masked_bases = int(line_parts[1])
            except ValueError:
                logger.error('line {}: the incorrect numeric value {'
                             '}'.format(self.__lineno, line_parts[1]))
                raise LavAlignmentError
        else:
            logger.error('line {}: an incorrect score record in an'
                         'x-stanza'.format(self.__lineno))

        # make sure that the stanza has the proper end
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        if self.__line != '}':
            logger.error('line {}: the x-stanza missed the proper '
                         'ending'.format(self.__lineno))
            raise LavAlignmentError

        return newly_masked_bases

    def __parse_m_stanza(self):
        """
        Parse the m-stanza from a LAV file.

        :return: the number of newly masked bases from the m-stanza
        :rtype: int
        """
        if self.__line != 'm {':
            logger.error('line {}: an incorrect m-stanza header')
            raise LavAlignmentError

        # the lines containing masked regions
        self.__line = self.__handler.readline().rstrip()
        self.__lineno += 1
        n_flag = False
        result_regions = []
        result_base_count = None
        while not n_flag:
            # we read masked region lines unless we read a line
            # starting with the n symbol
            line_parts = self.__line.split(None, 3)
            if line_parts[0] == 'm':
                if len(line_parts) == 3:
                    for i in xrange(1, 3):
                        try:
                            line_parts[i] = int(line_parts[i])
                        except ValueError:
                            logger.error(
                                'line {}: the incorrect numeric '
                                'value {}'.format(self.__lineno,
                                                  line_parts[i]))
                    result_regions.append(Lav.MStanzaItem(
                        start=line_parts[1], end=line_parts[2]))
                else:
                    logger.error('line {}: an incorrest masked '
                                 'region line'.format(self.__lineno))
                    raise LavAlignmentError
            elif line_parts[0] == 'n':
                if len(line_parts) == 2:
                    try:
                        result_base_count = int(line_parts[1])
                    except ValueError:
                        logger.error(
                            'line {}: the incorrect numeric value {'
                            '}'.format(self.__lineno, line_parts[1]))
                        raise ValueError
                    n_flag = True
                else:
                    logger.error(
                        'line {}: an incorrect total masked bases '
                        'line'.format(self.__lineno))
                    raise LavAlignmentError
            else:
                # the line inside an m stanza must start from either
                # 'x' or 'n', otherwise show the error message and
                # throw the exception
                logger.error('line {}: an incorrect line inside the m-'
                             'stanza'.format(self.__lineno))
                raise LavAlignmentError
            # read the next line
            self.__line = self.__handler.readline().rstrip()
            self.__lineno += 1

        result = Lav.MStanza(regions=result_regions,
                             base_count=result_base_count)

        # make sure that the stanza has the proper end
        if self.__line != '}':
            logger.error('line {}: the m-stanza missed the proper '
                         'ending'.format(self.__lineno))
            raise LavAlignmentError

        return result

    def __parse_census_stanza(self):
        """
        Parse a Census-stanza from a LAV file.
        
        :return: the list of Census-stanza contents
        :rtype: list
        """
        # check the first self.__line of the stanza
        if self.__line != 'Census {':
            logger.error('line {}: an incorrect Census-stanza '
                         'header'.format(self.__lineno))
            raise LavAlignmentError

        # start reading the Census stanza contents
        result = []
        i = 1
        while self.__line != '}':
            line_parts = self.__line.split(None, 2)
            try:
                line_parts = [int(x) for x in line_parts]
                if line_parts[0] != i:
                    logger.error(
                        'line {}: the incorrect Census record number {'
                        '}'.format(self.__lineno, i))
                    raise LavAlignmentError
                else:
                    result.append(line_parts[1])
            except ValueError:
                logger.error('line {}: the incorrect numeric value {'
                             '} or {}'.format(self.__lineno,
                                              *line_parts))
                raise LavAlignmentError
            i += 1
            # read the next self.__line
            self.__line = self.__handler.readline().rstrip()
            self.__lineno += 1

        return result
