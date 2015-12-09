|Travis| |PyPI| |Landscape| |Coveralls|

===========
Chromosomer
===========

Chromosomer is a reference-assisted assembly tool for producing draft
chromosome sequences. It provides the following routines:

 - ``fragmentmap`` - produce a fragment map from fragment alignments to reference chromosomes;

 - ``assemble`` - obtain FASTA sequences of assembled chromosomes from a fragment map;

 - ``transfer`` - move annotated regions from original fragments to assembled chromosomes;

 - ``fragmentmapstat`` - get summary on a fragment map;

 - ``fragmentmapbed`` - convert a fragment map to the BED format (e.g., for viewing in a genome browser);

 - ``fastalength`` - get lengths of sequences in a FASTA file (required for ``fragmentmap``).

Installation
------------

We recommend to install Chromosomer by using **pip**::

    pip install chromosomer

**pip** will automatically resolve Chromosomer dependencies and
install missing packages.

Assembling chromosomes
----------------------

A typical chromosome assembly involves two steps. First, a fragment
map is derived from fragment alignments to reference chromosome
sequences. Second, FASTA sequences of the assembled chromosomes are
obtained from the produced fragment map and the original fragment
sequences.

The first step is carried out with ``fragmentmap`` in the following
way::

    chromosomer fragmentmap alignment_file gap_size fragment_lengths output_map

The second step is implemented in the ``assemble`` routine::

    chromosomer assemble map fragment_fasta output_fasta

.. |PyPI| image:: https://img.shields.io/pypi/v/chromosomer.svg?branch=master
    :target: https://pypi.python.org/pypi/chromosomer
.. |Travis| image:: https://travis-ci.org/gtamazian/chromosomer.svg?branch=master
    :target: https://travis-ci.org/gtamazian/chromosomer
.. |Coveralls| image:: https://coveralls.io/repos/gtamazian/chromosomer/badge.svg?branch=master 
    :target: https://coveralls.io/r/gtamazian/chromosomer?branch=master
.. |Landscape| image:: https://landscape.io/github/gtamazian/chromosomer/master/landscape.svg?style=flat
   :target: https://landscape.io/github/gtamazian/chromosomer/master
   :alt: Code Health

