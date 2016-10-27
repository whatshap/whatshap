=======
Changes
=======

v0.13 (2016-10-27)
------------------

* Use ``PS`` tag instead of ``HP`` tag by default to store phasing information.
  This applies to the ``phase`` and ``hapcut2vcf`` subcommands. ``PS`` is also
  used by other tools and standard according to the VCF specification.
* Incorporated genotype likelihoods into our phasing framework. On request
  (by using option ``--distrust-genotypes``), genotypes can now be changed at a cost
  corresponding to their input genotype likelihoods. The changed genotypes are
  written to the output VCF. The behavior of ``--distrust-genotypes`` can be
  fine-tuned by the added options ``--include-homozygous``, ``--default-gq``,
  ``--gl-regularizer``, and ``--changed-genotype-list``.
* Correctly handle cases when processing VCFs with two or more disjoint
  families.

v0.12 (2016-07-01)
------------------

* Speed up allele detection
* Add an ``unphase`` subcommand which removes all phasing from a VCF file
  (``HP`` and ``PS`` tags, pipe notation).
* Add option ``--tag=`` to the ``phase`` subcommand, which allows to choose
  whether ReadBackedPhasing-compatible ``HP`` tags or standard ``PS`` tags are
  used to describe phasing in the output VCF.
* Manage versions with `versioneer <https://github.com/warner/python-versioneer>`_.
  This means that ``whatshap --version`` and the program version in the VCF header
  will include the Git commit hash, such as ``whatshap 0.11+50.g1b7af7a``.
* Add subcommand "haplotag" to tag reads in a BAM file with their haplotype.
* Fix a bug where re-alignment around variants at the very end of a chromosome
  would lead to an AssertionError.

v0.11 (2016-06-09)
------------------

* When phasing a pedigree, blocks that are not connected by reads but
  can be phased based on genotypes will be connected per default. This
  behavior can be turned off using option ``--no-genetic-haplotyping``.
* Implemented allele detection through re-alignment: To detect which allele of a
  variant is seen in a read, the query is aligned to the two haplotypes at that
  position. This results in better quality phasing, especially for
  low-quality reads (PacBio). Enabled if ``--reference`` is provided. Current
  limitation: No score for the allele is computed.
* As a side-effect of the new allele detection, we can now also phase
  insertions, deletions, MNPs and "complex" variants.
* Added option ``--chromosome`` to only work on specifed chromosomes.
* Use constant recombination rate per default, allows to use ``--ped``
  without using ``--genmap``.
* ``whatshap`` has become a command with subcommands. From now on, you need
  to run ``whatshap phase`` to phase VCFs.
* Add a ``stats`` subcommand that prints statistics about phased VCFs.

v0.10 (2016-04-27)
------------------

* Use ``--ped`` to phase pedigrees with the PedMEC algorithm
* Phase all samples in a multi-sample VCF
* Drop support for Python 3.2 - we require at least Python 3.3 now

v0.9 (2016-01-05)
-----------------

* This is the first release available via PyPI (and that can therefore be
  installed via ``pip install whatshap``)

January 2016
------------

* Trio phasing implemented in a branch

September 2015
--------------

* pWhatsHap implemented (in a branch)

April 2015
----------

* Create haplotype-specific BAM files

February 2015
-------------

* Smart read selection

January 2015
------------

* Ability to read multiple BAM files and merge them on the fly

December 2014
-------------

* Logo
* Unit tests

November 2014
-------------

* Cython wrapper for C++ code done
* Ability to write a phased VCF (using HP tags).

June 2014
---------

* Repository for WhatsHap refactoring created

April 2014
----------

* The WhatsHap algorithm is introduced at RECOMB
