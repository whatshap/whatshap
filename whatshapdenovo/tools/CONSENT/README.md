# CONSENT

CONSENT (Scalable long read self-correction and assembly polishing with multiple sequence alignment) is a self-correction method for long reads.
It works by, first, computing overlaps between the long reads, in order to define an alignment pile (i.e. a set of overlapping reads used for
correction) for each read. Each read's alignment pile is then further divided into smaller windows, that are corrected idependently.
First, a multiple alignment strategy is used in order to compute consensus. Then, this consensus is further polished with a local de Bruijn
graph, in order to get rid of the remaining errors.
Additionally to error correction, CONSENT can also perform assembly polishing.

Requirements
--------------

  - A Linux based operating system.
  - Python3.
  - g++, minimum version 5.5.0.
  - CMake, minimum version 2.8.2.
  - Minimap2 available through you path.
  
Installation
--------------

Clone the CONSENT repository, along with its submodules with:

  ```bash
  git clone --recursive https://github.com/morispi/CONSENT
  ```

Then run the install.sh script:

  ```bash
  ./install.sh
  ```

If you do not already have minimap2 available through your path, you can then run:

```bash
export PATH=$PWD/minimap2:$PATH
```

CONSENT should then be able to run.

Getting started
--------------

An example dataset (10x of simulated PacBio reads, raw assembly, and reference genome) is provided in the `example` folder.

Please run the following commands to try out CONSENT on this example.

### Self-correction

To perform self-correction on the example dataset, run the following command:

`./CONSENT-correct --in example/reads.fasta --out example/correctedReads.fasta --type PB`

This should take about 2 min and use up to 750 MB of RAM, using 4 cores.

### Polishing

To perform assembly polishing on the example dataset, run the following command:

`./CONSENT-polish --contigs example/rawAssembly.fasta --reads example/reads.fasta --out example/polishedAssembly.fasta`

This should take about 15 sec and use at most 150 MB of RAM, using 4 cores.


  
Running CONSENT
--------------

### Self-correction

To run CONSENT for long reads self-correction, run the following command:

`./CONSENT-correct --in longReads.fast[a|q] --out result.fasta --type readsTechnology`

  - longReads.fast[a|q]:	fasta or fastq file of long reads to .
  - result.fasta:		fasta file where to output the corrected long reads.
  - readsTechnology:	Indicate whether the long reads are from PacBio (--type PB) or Oxford Nanopore (--type ONT)


### Polishing

To run CONSENT for assembly polishing, run the followning command:

`./CONSENT-polish --contigs contigs.fast[a|q] --reads longReads.fast[a|q] --out result.fasta`

  - contigs.fast[a|q]:		fasta or fastq file of contigs to polish.
  - longReads.fast[a|q]:	fasta or fastq file of long reads to use for polishing.
  - result.fasta:		fasta file where to output the polished contigs.

### Options

      --windowSize INT, -l INT:      Size of the windows to process. (default: 500)
      --minSupport INT, -s INT:      Minimum support to consider a window for correction. (default: 4)
      --maxSupport INT, -S INT:      Maximum number of overlaps to include in a pile. (default: 150)
      --maxMSA INT, -M:              Maximum number of sequences to include into the MSA. (default: 150)
      --merSize INT, -k INT:         k-mer size for chaining and polishing. (default: 9)
      --solid INT, -f INT:           Minimum number of occurrences to consider a k-mer as solid during polishing. (default: 4)
      --anchorSupport INT, -c INT:   Minimum number of sequences supporting (Ai) - (Ai+1) to keep the two anchors in the chaining. (default: 8)
      --minAnchors INT, -a INT:      Minimum number of anchors in a window to allow consensus computation. (default: 2)
      --windowOverlap INT, -o INT:   Overlap size between consecutive windows. (default: 50)
      --nproc INT, -j INT:           Number of processes to run in parallel (default: number of cores).
      --minimapIndex INT, -m INT:    Split minimap2 index every INT input bases (default: 500M).
      --tmpdir STRING, -t STRING:    Path where to store the temporary overlaps file (default: working directory, as Alignments_dateTimeStamp.paf).
      --help, -h:                    Print this help message.

Notes
--------------

CONSENT has been developed and tested on x86-64 GNU/Linux.          
Support for any other platform has not been tested.

Authors
--------------

Pierre Morisse, Camille Marchet, Antoine Limasset, Arnaud Lefebvre and Thierry Lecroq.

Reference
--------------

Morisse, P., Marchet, C., Limasset, A. et al. Scalable long read self-correction and assembly polishing with multiple sequence alignment. Sci Rep 11, 761 (2021). https://doi.org/10.1038/s41598-020-80757-5

Contact
--------------

You can report problems and bugs to pierre[dot]morisse[at]inria[dot]fr
