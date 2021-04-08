.. image:: https://img.shields.io/pypi/v/whatshap.svg?branch=main
    :target: https://pypi.python.org/pypi/whatshap

.. image:: https://github.com/whatshap/whatshap/workflows/CI/badge.svg

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
    :target: http://bioconda.github.io/recipes/whatshap/README.html

|

.. image:: https://github.com/whatshap/whatshap/raw/main/logo/whatshap_logo.png

|

WhatsHap
========

WhatsHap is a software for phasing genomic variants using DNA sequencing
reads, also called *read-based phasing* or *haplotype assembly*. It is
especially suitable for long reads, but works also well with short reads.


Features
========

  * Very accurate results (Martin et al.,
    `WhatsHap: fast and accurate read-based phasing <https://doi.org/10.1101/085050>`_)
  * Works well with Illumina, PacBio, Oxford Nanopore and other types of reads
  * It phases SNVs, indels and even “complex” variants (such as ``TCG`` → ``AGAA``)
  * Pedigree phasing mode uses reads from related individuals (such as trios)
    to improve results and to reduce coverage requirements
    (Garg et al., `Read-Based Phasing of Related Individuals <https://doi.org/10.1093/bioinformatics/btw276>`_).
  * WhatsHap is easy to install
  * It is easy to use: Pass in a VCF and one or more BAM files, get out a phased VCF.
    Supports multi-sample VCFs.
  * It produces standard-compliant VCF output by default
  * If desired, get output that is compatible with ReadBackedPhasing
  * Open Source (MIT license)


Documentation
-------------

* `GitHub repository <https://github.com/whatshap/whatshap/>`_
* `Read the documentation online <https://whatshap.readthedocs.io/>`_.
  Offline documentation is available in the ``doc/`` subdirectory in the
  repository and in the downloaded tar distribution.


Issue tracker
-------------
Please do not hesitate to use our `issue tracker <https://github.com/whatshap/whatshap/issues>`_ for bug reports and feature requests.
