#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com


class AlignmentError(Exception):
    """
    The class that describes a basic error that may occur while
    parsing an alignment file.
    """
    pass


class BlastAlignmentError(AlignmentError):
    """
    The class describes an error that may occur while parsing a file
    in the BLAST tabular output format.
    """
    pass


class LavAlignmentError(AlignmentError):
    """
    The class describes an error that may occur while parsing a file
    in the LAV format.
    """
    pass
