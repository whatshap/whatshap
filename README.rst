.. image:: https://img.shields.io/pypi/v/whatshap.svg?branch=master
    :target: https://pypi.python.org/pypi/whatshap

.. image:: https://semaphoreci.com/api/v1/whatshap/whatshap/branches/master/shields_badge.svg
    :target: https://semaphoreci.com/whatshap/whatshap

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
    :target: http://bioconda.github.io/recipes/whatshap/README.html

|

.. image:: https://bitbucket.org/repo/8AjxBd/images/543323972-whatshap_logo.png
    :scale: 50%

|

WhatsHap
========

WhatsHap is a software for phasing genomic variants using DNA sequencing
reads, also called *read-based phasing* or *haplotype assembly*. It is
especially suitable for long reads, but works also well with short reads.

Please :ref:`cite us if you use WhatsHap <howtocite>`.


Features
========

  * WhatsHap is :ref:`easy to install <installation>`
  * It is :ref:`easy to use <user-guide>`: Pass in a VCF and a BAM, get out a phased VCF
  * Works well with Illumina, PacBio, Oxford Nanopore and other types of reads
  * It phases indels
  * It produces standard-compliant VCF output (or optionally output that is compatible with ReadBackedPhasing)
  * Highly accurate results (Martin et al.,
    `WhatsHap: fast and accurate read-based phasing <https://doi.org/10.1101/085050>`_)
  * Phasing of pedigrees: If you have reads from a trio or any pedigree of
    related individuals, phase them at the same time in pedigree mode to improve
    results and to lower coverage requirements at the same time
    (Garg et al., `Read-Based Phasing of Related Individuals <https://doi.org/10.1093/bioinformatics/btw276>`_).
  * Open Source (MIT license)


Documentation
-------------

* `Bitbucket page <https://bitbucket.org/whatshap/whatshap/>`_
* `Read the documentation online <https://whatshap.readthedocs.io/>`_.
  Offline documentation is available in the ``doc/`` subdirectory in the
  repository and in the downloaded tar distribution.


Mailing list
------------
We run a `public mailing list <https://lists.cwi.nl/mailman/listinfo/whatshap>`_. Please
don't hesitate to post questions and comments.
