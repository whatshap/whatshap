.. image:: https://img.shields.io/pypi/v/whatshap.svg?branch=master
    :target: https://pypi.python.org/pypi/whatshap

.. image:: https://semaphoreci.com/api/v1/whatshap/whatshap/branches/master/shields_badge.svg
    :target: https://semaphoreci.com/whatshap/whatshap

|

.. image:: https://bitbucket.org/repo/8AjxBd/images/3378940113-whatshap_logo.png
    :height: 150px

WhatsHap
========

WhatsHap is a software for phasing genomic variants using DNA sequencing
reads, also called *haplotype assembly*. It is especially suitable for long
reads, but works also well with short reads.

If you use WhatsHap, please cite:

    Murray Patterson, Tobias Marschall, Nadia Pisanti, Leo van Iersel,
    Leen Stougie, Gunnar W. Klau, Alexander Schönhuth.
    `WhatsHap: Weighted Haplotype Assembly for Future-Generation Sequencing Reads
    <http://dx.doi.org/10.1089/cmb.2014.0157>`_.
    Journal of Computational Biology, 22(6), pp. 498-509, 2015.
    (`Get a self-archived version <https://bioinf.mpi-inf.mpg.de/homepage/publications.php?&account=marschal>`_)

The version of WhatsHap you find here is the result of further development
focused on making the software easy and straightforward to use. WhatsHap is now
Open Source software under the MIT license and we welcome contributions.


Pedigree phasing
----------------

WhatsHap is capable of :ref:`using pedigree information <phasing-pedigrees>`
about individuals to further improve phasing results, and to drastically reduce
the needed coverage. A preprint is available on bioRxiv:

    Read-Based Phasing of Related Individuals.
    Shilpa Garg, Marcel Martin, Tobias Marschall.
    `doi:10.1101/037101 <http://dx.doi.org/10.1101/037101>`_


Parallel Version: pWhatsHap
---------------------------
A parallelization of the core dynamic programming algorithm has been described in 

    M. Aldinucci, A. Bracciali, T. Marschall, M. Patterson, N. Pisanti, M. Torquati. 
    `High-Performance Haplotype Assembly <http://dx.doi.org/10.1007/978-3-319-24462-4_21>`_. Proceedings of the 11th International
    Meeting on Computational Intelligence Methods for Bioinformatics and
    Biostatistics (CIBB), 245-258, 2015.

The current implementation can be found in `branch parallel <https://bitbucket.org/whatshap/whatshap/branch/parallel>`_.


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
