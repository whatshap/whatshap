Parts of WhatsHap have been described in different articles. Please choose
an appropriate citation depending on your use case.

If you use WhatsHap as a tool:

    | Marcel Martin, Murray Patterson, Shilpa Garg, Sarah O. Fischer,
      Nadia Pisanti, Gunnar W. Klau, Alexander Schoenhuth, Tobias Marschall.
    | *WhatsHap: fast and accurate read-based phasing*
    | bioRxiv 085050
    | doi: `10.1101/085050 <https://doi.org/10.1101/085050>`_

To refer to the core WhatsHap phasing algorithm:

    | Murray Patterson, Tobias Marschall, Nadia Pisanti, Leo van Iersel,
      Leen Stougie, Gunnar W. Klau, Alexander Schönhuth.
    | *WhatsHap: Weighted Haplotype Assembly for Future-Generation Sequencing Reads*
    | Journal of Computational Biology, 22(6), pp. 498-509, 2015.
    | doi: `10.1089/cmb.2014.0157 <http://dx.doi.org/10.1089/cmb.2014.0157>`_
      (`Get self-archived PDF <https://bioinf.mpi-inf.mpg.de/homepage/publications.php?&account=marschal>`_)

To refer to the pedigree-phasing algorithm and the PedMEC problem:

    | Shilpa Garg, Marcel Martin, Tobias Marschall.
    | *Read-based phasing of related individuals*
    | Bioinformatics 2016; 32 (12): i234-i242.
    | doi: `10.1093/bioinformatics/btw276 <https://doi.org/10.1093/bioinformatics/btw276>`_

WhatsHap's genotyping algorithm is described here:

    | Jana Ebler, Marina Haukness, Trevor Pesout, Tobias Marschall, Benedict Paten.
    | *Haplotype-aware genotyping from noisy long reads*
    | bioRxiv
    | doi: `10.1101/293944 <https://doi.org/10.1101/293944>`_

The HapChat algorithm is an alternative MEC solver able to handle higher coverages. It can be used
through "whatshap phase --algorithm=hapchat". It has been described in this paper:

    | Stefano Beretta, Murray Patterson, Simone Zaccaria, Gianluca Della Vedova, Paola Bonizzoni.
    | *HapCHAT: adaptive haplotype assembly for efficiently leveraging high coverage in long reads*.
    | BMC Bioinformatics, 19:252, 2018.
    | doi: `10.1186/s12859-018-2253-8 <https://doi.org/10.1186/s12859-018-2253-8>`_
    
A parallelization of the core dynamic programming algorithm (“pWhatsHap”)
has been described in

    | M. Aldinucci, A. Bracciali, T. Marschall, M. Patterson, N. Pisanti, M. Torquati.
    | *High-Performance Haplotype Assembly*
    | Proceedings of the 11th International Meeting on Computational Intelligence
      Methods for Bioinformatics and Biostatistics (CIBB), 245-258, 2015.
    | doi: `10.1007/978-3-319-24462-4_21 <http://dx.doi.org/10.1007/978-3-319-24462-4_21>`_

pWhatsHap is currently not integrated into the main WhatsHap source code. It
is available in
`branch parallel <https://bitbucket.org/whatshap/whatshap/branch/parallel>`_
in the Git repository.

If you use the polyploid phasing algorithm (``whatshap polyphase``), please refer to

    | Sven D. Schrinner, Rebecca Serra Mari, Jana Ebler, Mikko Rautiainen, Lancelot Seillier,
    | Julia J. Reimer, Björn Usadel, Tobias Marschall, Gunnar W. Klau.
    | *Haplotype threading: accurate polyploid phasing from long reads*
    | Genome Biology
    | doi: `10.1186/s13059-020-02158-1 <https://doi.org/10.1186/s13059-020-02158-1>`_

If you use the polyploid phasing algorithm that takes pedigree information into account
(``whatshap geneticpolyphase``), plese refer to

    | Sven Schrinner, Rebecca Serra Mari, Richard Finkers, Paul Arens,
    | Björn Usadel, Tobias Marschall, Gunnar W. Klau.
    | *Genetic polyploid phasing from low-depth progeny samples*.
    | iScience 25(6), 2022.
    | doi: `10.1016/j.isci.2022.104461 <https://doi.org/10.1016/j.isci.2022.104461>`_
