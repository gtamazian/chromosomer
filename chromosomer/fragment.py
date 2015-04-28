#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import pyfaidx
import random
import string
from chromosomer.alignment.blast import Blast
from chromosomer.exception import MapError
from chromosomer.exception import AlignmentToMapError
from chromosomer.fasta import RandomSequence
from chromosomer.fasta import Writer
from collections import defaultdict
from collections import namedtuple
from operator import attrgetter

logging.basicConfig()
logger = logging.getLogger(__name__)


class Map(object):
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
        Initializes a Map object.
        """
        self.__fragments = defaultdict(list)
        self.__block_adding = False

    def add_record(self, new_record):
        """
        Given a new fragment record, add it to the fragment map.

        :param new_record: a record to be added to the map
        :type new_record: Map.Record
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
                    raise MapError
                for i in self.numeric_values:
                    try:
                        line_parts[i] = int(line_parts[i])
                    except ValueError:
                        logger.error('line %d: the incorrect numeric '
                                     'value %s', lineno, line_parts[i])
                        raise MapError
                new_record = Map.Record(*line_parts)
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
            raise MapError

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
        with Writer(output_filename) as chromosome_writer:
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
                            raise MapError
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


class AlignmentToMap(object):
    """
    The class implements routines to create a fragment map from a set
    of alignments between fragments to be assembled and reference
    chromosomes.
    """

    Anchor = namedtuple('Anchor', ('fragment', 'fr_start', 'fr_end',
                                   'fr_strand', 'ref_chr',
                                   'ref_start', 'ref_end'))

    def __init__(self, gap_size, fragment_lengths,
                 min_fragment_length=None, centromeres=None):
        """
        Create a converter object to create fragment maps from
        alignmets between reference chromosomes and fragments to be
        assembled.

        :param gap_size: a size of a gap between fragments
        :param fragment_lengths: a dictionary of fragment lengths
            which keys are their names and values are their lengths
        :param min_fragment_length: the minimal length of a fragment
            to be included in the map
        :param centromeres: a dictionary of reference chromosome
            centromere locations
        :type gap_size: int
        :type fragment_lengths: dict
        :type min_fragment_length: int
        :type centromeres: dict
        """
        self.__gap_size = gap_size
        self.__fragment_lengths = fragment_lengths
        self.__min_fragment_length = min_fragment_length
        self.__centromeres = centromeres

        self.__anchors = {}
        self.__unlocalized = []
        self.__unplaced = []
        self.__fragment_map = Map()

    def blast(self, blast_alignments, bitscore_ratio_threshold):
        """
        Create a fragment map from BLAST blast_alignments between
        fragments and reference chromosomes.

        :param blast_alignments: BLAST blast_alignments
        :param bitscore_ratio_threshold:
        :type blast_alignments: Blast
        :return: the fragment map constructed from the blast_alignments
        :rtype: Map
        """
        self.__anchors = {}
        self.__unlocalized = []
        self.__unplaced = []

        temp_anchors = defaultdict(list)

        for alignment in blast_alignments.alignments():
            if self.__min_fragment_length is not None:
                # check if the fragment length is equal or greater
                # than the threshold value
                try:
                    if self.__fragment_lengths[alignment.query] < \
                            self.__min_fragment_length:
                        # skip the alignment
                        continue
                except KeyError:
                    logger.error('the fragment %s length is missing',
                                 alignment.query)
                    raise AlignmentToMapError

            # consider the centromeres if required
            if self.__centromeres is not None and alignment.subject \
                    in self.__centromeres:
                # the chromosome a fragment was aligned to has a
                # centromere, so we determine which arm the alignment
                # refers to and modify the chromosome name by adding
                # '_1' or '_2' to it
                if min(alignment.s_start, alignment.s_end) < \
                        self.__centromeres[alignment.subject].start:
                    arm_prefix = '_1'
                else:
                    arm_prefix = '_2'
                new_alignment = list(alignment)
                new_alignment[1] += arm_prefix
                alignment = Blast.Alignment(*new_alignment)

            temp_anchors[alignment.query].append(alignment)
            # check if there is more than 2 alignments for the
            # fragment; if there is, then leave two fragments with
            # the greatest bit-score values
            if len(temp_anchors[alignment.query]) > 2:
                temp_anchors[alignment.query] = sorted(
                    temp_anchors[alignment.query],
                    key=attrgetter('bit_score'),
                    reverse=True
                )
                temp_anchors[alignment.query] = temp_anchors[
                    alignment.query][0:2]

        for fragment, alignments in temp_anchors.iteritems():
            if len(alignments) > 1:
                # check if the ratio of the alignment bit scores is
                # greater than the required threshold to consider a
                # fragment places
                if alignments[0].bit_score/alignments[1].bit_score > \
                        bitscore_ratio_threshold:
                    self.__anchors[fragment] = \
                        AlignmentToMap.Anchor(
                            fragment=alignments[0].fragment,
                            fr_start=alignments[0].fr_start - 1,
                            fr_end=alignments[0].fr_end,
                            fr_strand=alignments[0].fr_strand,
                            ref_chr=alignments[0].ref_chr,
                            ref_start=min(alignments[0].ref_start,
                                          alignments[0].ref_end) - 1,
                            ref_end=max(alignments[0].ref_start,
                                        alignments[0].ref_end)
                        )
                elif alignments[0].subject == alignments[1].subject:
                    # the fragment is considered unlocalized
                    self.__unlocalized.append(
                        (fragment, alignments[0].subject))
                else:
                    # the fragment is considered unplaced
                    self.__unplaced.append(fragment)
            else:
                # there is a single alignment, use it as an anchor
                self.__anchors[fragment] = \
                    AlignmentToMap.Anchor(
                        fragment=alignments[0].fragment,
                        fr_start=alignments[0].fr_start - 1,
                        fr_end=alignments[0].fr_end,
                        fr_strand=alignments[0].fr_strand,
                        ref_chr=alignments[0].ref_chr,
                        ref_start=min(alignments[0].ref_start,
                                      alignments[0].ref_end) - 1,
                        ref_end=max(alignments[0].ref_start,
                                    alignments[0].ref_end)
                    )

        self.__anchor_fragments()

    def __anchor_fragments(self):
        """
        Build a fragment map from anchors.
        """
        # first, we split anchors by reference genome chromosomes
        chr_anchors = defaultdict(list)
        for anchor in self.__anchors.itervalues():
            chr_anchors[anchor.subject].append(anchor)

        # second, we sort the anchors by their position on the
        # chromosomes
        for chr_name in chr_anchors.iterkeys():
            chr_anchors[chr_name] = sorted(
                chr_anchors[chr_name],
                key=lambda x: min(x.s_start, x.s_end)
            )

        # now we form a fragment map from the anchors
        self.__fragment_map = Map()
        for chr_name in chr_anchors.iterkeys():
            for anchor in chr_anchors[chr_name]:
                try:
                    fragment_length = self.__anchors[anchor.query]
                except ValueError:
                    logger.error('the fragment %s length is missing',
                                 anchor.query)
                    raise AlignmentToMapError

                # determine the fragment's start and end positions
                if anchor.s_start < anchor.s_end:
                    fragment_strand = '+'
                    ref_start = anchor.ref_start - anchor.fr_start
                    ref_end = ref_start + fragment_length
                else:
                    fragment_strand = '-'
                    ref_end = anchor.ref_end + anchor.fr_start
                    ref_start = ref_end - fragment_length

                new_record = Map.Record(
                    fr_name=anchor.query,
                    fr_length=fragment_length,
                    fr_start=0,
                    fr_end=fragment_length,
                    fr_strand=fragment_strand,
                    ref_chr=anchor.subject,
                    ref_start=ref_start,
                    ref_end=ref_end
                )
                self.__fragment_map.add_record(new_record)


class Simulator(object):
    """
    The class describes routines to simulate genome fragments and
    chromosomes that are composed from them.
    """
    def __init__(self, fragment_length, fragment_number,
                 chromosome_number, gap_size):
        """
        Create a fragment simulator object.

        :param fragment_length: the length of a fragment
        :param fragment_number: the number of fragments constituting
            the chromosomes
        :param chromosome_number: the number of chromosomes
        :param gap_size: the length of gaps between fragments in
            chromosomes
        :type fragment_length: int
        :type fragment_number: int
        :type chromosome_number: int
        :type gap_size: int
        """
        self.__fragment_length = fragment_length
        self.__fragment_number = fragment_number
        self.__chromosome_number = chromosome_number
        self.__gap_size = gap_size

        # create fragment sequences
        self.__fragments = {}
        seq_generator = RandomSequence(self.__fragment_length)
        for i in xrange(self.__fragment_number):
            fr_name = 'fragment{}'.format(i+1)
            self.__fragments[fr_name] = seq_generator.get()

        self.__map = Map()
        self.__create_map()
        self.__assemble_chromosomes()

    def __create_map(self):
        """
        Assign fragments to chromosomes randomly and create a
        fragment map.
        """
        fragment_positions = [0] * self.__chromosome_number
        for i in xrange(self.__fragment_number):
            chr_num = random.randrange(self.__chromosome_number)
            fr_strand = random.choice(('+', '-'))
            self.__map.add_record(Map.Record(
                fr_name='fragment{}'.format(i+1),
                fr_length=self.__fragment_length,
                fr_start=0,
                fr_end=self.__fragment_length,
                fr_strand=fr_strand,
                ref_chr='chr{}'.format(chr_num+1),
                ref_start=fragment_positions[chr_num],
                ref_end=fragment_positions[chr_num] +
                self.__fragment_length
            ))
            fragment_positions[chr_num] += self.__gap_size
            self.__map.add_record(Map.Record(
                fr_name='GAP',
                fr_length=self.__gap_size,
                fr_start=0,
                fr_end=self.__gap_size,
                fr_strand='+',
                ref_chr='chr{}'.format(chr_num+1),
                ref_start=fragment_positions[chr_num],
                ref_end=fragment_positions[chr_num] + self.__gap_size
            ))
            fragment_positions[chr_num] += self.__gap_size

    def __assemble_chromosomes(self):
        """
        Get chromosome sequences from fragments using the constructed
        fragment map.
        """
        complement = string.maketrans('ATCGatcgNnXx', 'TAGCtagcNnXx')
        chromosomes = defaultdict(list)
        for i in self.__map.chromosomes():
            for fr in self.__map.fragments(i):
                if fr.fr_name == 'GAP':
                    temp_fragment = 'N' * fr.fr_length
                else:
                    temp_fragment = self.__fragments[fr.fr_name]
                    if fr.fr_strand == '-':
                        temp_fragment = temp_fragment.translate(
                            complement)
                chromosomes[i].append(temp_fragment)
            chromosomes[i] = ''.join(chromosomes[i])

        self.__chromosomes = chromosomes

    def write(self, map_file, fragment_file, chromosome_file):
        """
        Write the produced data - a fragment map, a FASTA file of
        fragments and a FASTA file of chromosomes - to the specified
        files.

        :param map_file: a name of a file to write the fragment map to
        :param fragment_file: a name of a file to write fragment
            sequences to
        :param chromosome_file: a name of a file to write chromosome
            sequences to
        :type map_file: str
        :type fragment_file: str
        :type chromosome_file: str
        """
        self.__map.write(map_file)
        with Writer(fragment_file) as fragment_fasta:
            for i, seq in enumerate(self.__fragments):
                fragment_fasta.write('fragment{}'.format(i+1), seq)
        with Writer(chromosome_file) as chromosome_fasta:
            for i, seq in enumerate(self.__chromosomes):
                chromosome_fasta.write('chr{}'.format(i+1), seq)


class Length(object):
    """
    The class implements routines to handle fragment sequence lengths.
    """
    def __init__(self, filename):
        """
        Create a Length object to handle sequence lengths of the
        specified FASTA file.

        :param filename: a name of a FASTA file with sequences which
            lengths are to be derived
        :type filename: str
        """
        self.__filename = filename
        self.__lengths = {}

    def lengths(self):
        """
        Return a dictionary of sequence lengths.

        :return: a dictionary which keys are sequence names and
            values are their lengths
        :rtype: dict
        """
        if not self.__lengths:
            reader = pyfaidx.Fasta(self.__filename)
            for seq in reader.keys():
                self.__lengths[seq] = len(reader[seq])

        return self.__lengths
