#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import random


class Writer(object):
    """
    The class implements routines to write sequences in the FASTA
    format.
    """

    def __init__(self, filename, width=72):
        """
        Create a Writer object to write sequences in a FASTA file.

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


class RandomSequence(object):
    """
    The class implements routines to create random nucleotide
    sequences.
    """

    def __init__(self, length):
        """
        Create an object to generate random nucleotide sequences of
        the specified length.

        :param length: a sequence length
        :type length: int
        """
        self.__length = length

    def get(self):
        """
        Get a random nucleotide sequence.

        :return: a random nucleotide sequence
        :rtype: str
        """
        result = []
        nucleotides = ('A', 'C', 'G', 'T')
        for _ in xrange(self.__length):
            result.append(random.choice(nucleotides))

        return ''.join(result)
