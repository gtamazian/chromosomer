#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import logging
import subprocess
from itertools import chain

logging.basicConfig()
logger = logging.getLogger(__name__)


class BlastN(object):
    """
    The class implements a wrapper to launch blastn from the NCBI
    BLAST+ package.
    """

    def __init__(self, query, database, output):
        """
        Create a BlastN object to align the specified query to the
        specified database.

        :param query: a name of a FASTA file of query sequences to be
            aligned
        :param database: a name of a BLAST database to align the query
            sequences to
        :type query: str
        :type database: str
        """
        self.__query = query
        self.__database = database
        self.__output = output
        self.__parameters = {}

    def get(self, parameter):
        """
        Get a value of the specified parameter.

        :param parameter: a parameter name
        :type parameter: str
        :return: the specified parameter value or None if it was not
            specified
        """
        return self.__parameters.setdefault(parameter)

    def set(self, parameter, value):
        """
        Set the value of a blastn option.

        :param parameter: a parameter name
        :param value: a parameter value
        :type parameter: str
        """
        self.__parameters[parameter] = value

    def launch(self):
        """
        Launch blastn with the specified parameters.
        """
        options = ['blastn', '-query', self.__query, '-db',
                   self.__database, '-out', self.__output] + list(
            chain.from_iterable(self.__parameters.iteritems()))

        subprocess.check_call(options)
