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
    `WhatsHap: Haplotype Assembly for Future-Generation Sequencing
    Reads <http://dx.doi.org/10.1007/978-3-319-05269-4_19>`_.
    Proceedings of ACM 18th Annual International Conference on Research in
    Computational Biology (RECOMB), 237-249, 2014. 
    (A self-archived version can be found `here <https://bioinf.mpi-inf.mpg.de/homepage/publications.php?&account=marschal>`_)

The version of WhatsHap you find here is the result of further development
focused on making the software easy and straightforward to use. WhatsHap is now
Open Source software under the MIT license and we welcome contributions.


.. note:: WhatsHap is work in progress! In particular, the documentation is
	incomplete, not all features that we would like to have for an initial
	release are there, and there are probably bugs.

Parallel Version: pWhatsHap
---------------------------
A parallelization of the core dynamic programming algorithm has been described in 

    M. Aldinucci, A. Bracciali, T. Marschall, M. Patterson, N. Pisanti, M. Torquati. 
    `High-Performance Haplotype Assembly <http://dx.doi.org/10.1007/978-3-319-24462-4_21>`_. Proceedings of the 11th International
    Meeting on Computational Intelligence Methods for Bioinformatics and
    Biostatistics (CIBB), 245-258, 2015.

The current implementation can be found in `branch parallel <https://bitbucket.org/whatshap/whatshap/branch/parallel>`_.
   
Links
-----

* `Bitbucket page <https://bitbucket.org/whatshap/whatshap/>`_
* `Read the documentation online <https://whatshap.readthedocs.org/>`_.
  Offline documentation is available in the ``doc/`` subdirectory in the
  repository and in the downloaded tar distribution.