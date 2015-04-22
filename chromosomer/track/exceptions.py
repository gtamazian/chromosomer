#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com


class TrackError(Exception):
    """
    The class describes a basic error that may occur while reading or
    writing a track file.
    """
    pass


class BedError(TrackError):
    """
    The class describes an error that may occur while reading or
    writing a track file in the BED format.
    """
    pass
