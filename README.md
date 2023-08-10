[![PyPI](https://img.shields.io/pypi/v/whatshap.svg)](https://pypi.python.org/pypi/whatshap)
![CI](https://github.com/whatshap/whatshap/workflows/CI/badge.svg)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/whatshap/README.html)

![WhatsHap logo](https://github.com/whatshap/whatshap/raw/main/logo/whatshap_logo.png)

WhatsHap is a software for phasing genomic variants using DNA sequencing
reads, also called *read-based phasing* or *haplotype assembly*. It is
especially suitable for long reads, but works also well with short reads.

For documentation and information on how to cite WhatsHap, please visit the [WhatsHap Homepage](https://whatshap.readthedocs.io/)
# _k_-merald
This branch implements k-merald, an allele detection approach using k-mer based sequencing error profiles, as an alternative to the edit distance based allele detection implemented in the original WhatsHap version.
For installation, pelase follow the WhatsHap development version installation istructions available at [WhatsHap Homepage](https://whatshap.readthedocs.io/)
## Usage
**Note**: All the steps mentioned below should be run from within the conda environment created during the installation.

1) Train the model on non-variant regions of the genome:
   
   i) Get reference-read kmer pair counts using the **call** module
     ```bash
     whatshap call <reference_genome.fa> <aligned_reads.bam> <variants.vcf> <kmer_size> <window_size/2> <pseudocount_value_for_unobserved_kmer_pairs> > <ref-read_kmer_pair_counts>
     ```
     For window_size, you need to specify the window length you require on each side of the variant i.e. half of your required total window size. 
     **Note**: It is suggested to perform this training step in parallel on chromosome specific reads to save time, however, it is not mandatory.

   ii) Generate phred-scores using the kmer-pair counts
     ```bash
     python3 whatshap/phred_scores.py -i <ref-read_kmer_counts_dir> -o <phred_scores> -k <kmer_size> -e <pseudocount_value_for_unobserved_kmer_pairs>
     ```
     **Note**: This script takes the path of the folder containg all the output files from multiple chromosome specific iterations of step (i).

2) For genotyping use the **genotype** module 
    ```bash
     whatshap genotype [options] --reference <reference_genome.fa> -o <genotyped_variants.vcf> <variant_to_genotype.vcf> <aligned_reads.bam> -b <phred_scores> -k <kmer_size> -g {gap_value_for_kmer_alignment}
    ```
    
2) For phasing use the **phase** module 
    ```bash
     whatshap phase [options] --reference <reference_genome.fa> -o <genotyped_variants.vcf> <variant_to_genotype.vcf> <aligned_reads.bam> -b <phred_scores> -k <kmer_size> -g {gap_value_for_kmer_alignment}
    ```
    

   


   
     
