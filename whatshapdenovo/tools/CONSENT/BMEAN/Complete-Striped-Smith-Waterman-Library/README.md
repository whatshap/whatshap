# SSW Library

## An SIMD Smith-Waterman C/C++/Python/Java Library for Use in Genomic Applications

License: MIT

Copyright (c) 2012-2015 Boston College

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Authors:
* C implementation: Mengyao Zhao
* C++ wrapper: Wan-Ping Lee
* Python wrapper: Yongan Zhao
* Java wrapper: Daniel Cameron

Contact:
* Mengyao Zhao <zhaomengyao@gmail.com>
* Wan-Ping Lee <wanping.lee@gmail.com>
* Yongan Zhao <zhaoyanswill@gmail.com>
* Daniel Cameron <cameron.d@wehi.edu.au>

Last revision: 06/29/2016

## Overview

SSW is a fast implementation of the Smith-Waterman algorithm, which uses the Single-Instruction Multiple-Data (SIMD) instructions to parallelize the algorithm at the instruction level. SSW library provides an API that can be flexibly used by programs written in C, C++ and other languages. We also provide a software that can do protein and genome alignment directly. Current version of our implementation is ~50 times faster than an ordinary Smith-Waterman. It can return the Smith-Waterman score, alignment location and traceback path (cigar) of the optimal alignment accurately; and return the sub-optimal alignment score and location heuristically.

The Debian package of this library can be achieved here:
https://tracker.debian.org/pkg/libssw

Note: When SSW open a gap, the gap open penalty alone is applied.

## C/C++ interface

### How to use the API

The API files include `ssw.h` and `ssw.c`, which can be directly used by any C or C++ program. For the C++ users who are more comfortable to use a C++ style interface, an additional C++ wrapper is provided with the file `ssw_cpp.cpp` and `ssw_cpp.h`.

To use the C style API, please:

1. Download `ssw.h` and `ssw.c`, and put them in the same folder of your own program files.
2. Write `#include "ssw.h"` into your file that will call the API functions.
3. The API files are ready to be compiled together with your own C/C++ files.

The API function descriptions are in the file `ssw.h`. One simple example of the API usage is `example.c`. The Smith-Waterman penalties need to be integers. Small penalty numbers such as: match: 2, mismatch: -1, gap open (the total penalty when one gap is opened): -3, gap extension: -1 are recommended, which will lead to shorter running time.

To use the C++ style API, please:

1. Download `ssw.h`, `ssw.c`, `ssw_cpp.cpp` and `ssw_cpp.h`, put them in the same folder of your own program files.
2. Write `#include "ssw_cpp.h"` into your file that will call the API functions.
3. The API files are ready to be compiled together with your own C/C++ files.

The API function descriptions are in the file `ssw_cpp.h`. A simple example of using the C++ API is `example.cpp`.

### Speed and memory usage of the API

Test data set:

* Target sequence: reference genome of _E. coli_ strain 536 (4,938,920 nucleotides) from NCBI
* Query sequences: 1000 reads of Ion Torrent sequenced _E. coli_ strain DH10B (C23-140, 318 PGM Run, 11/2011), read length: ~25-540 bp, most reads are ~200 bp

CPU time:

* AMD CPU: default penalties: ~880 seconds; -m1 -x3 -o5 -e2: ~460 seconds
* Intel CPU: default penalties: ~960 seconds; -m1 -x3 -o5 -e2: ~500 seconds

Memory usage: ~40MB

### Install the software

1. Download the software from https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.
2. `cd src`
3. `make`
4. The executable file will be `ssw_test`.

### Run the software

```
Usage: ssw_test [options] ... <target.fasta> <query.fasta>(or <query.fastq>)
Options:
	-m N	N is a positive integer for weight match in genome sequence alignment. [default: 2]
	-x N	N is a positive integer. -N will be used as weight mismatch in genome sequence alignment. [default: 2]
	-o N	N is a positive integer. -N will be used as the weight for the gap opening. [default: 3]
	-e N	N is a positive integer. -N will be used as the weight for the gap extension. [default: 1]
	-p	Do protein sequence alignment. Without this option, the ssw_test will do genome sequence alignment.
	-a FILE	FILE is either the Blosum or Pam weight matrix. [default: Blosum50]
	-c	Return the alignment path.
	-f N	N is a positive integer. Only output the alignments with the Smith-Waterman score >= N.
	-r	The best alignment will be picked between the original read alignment and the reverse complement read alignment.
	-s	Output in SAM format. [default: no header]
	-h	If -s is used, include header in SAM output.
```

### Software input

The input files can be in FASTA or FASTQ format. Both target and query files can contain multiple sequences. Each sequence in the query file will be aligned with all sequences in the target file. If your target file has N sequences and your query file has M sequences, the results will have M*N alignments.

### Software output

The software can output SAM format or BLAST like format results.

1. SAM format output:

Example:

```
@HD VN:1.4  SO:queryname
@SQ SN:chr1 LN:1001
6:163296599:F:198;None;None/1   0   chr1    453 5   3M2D3M1D4M2D6M1D5M1D5M2I7M  *   0   0   CCAGCCCAAAATCTGTTTTAATGGTGGATTTGTGT *   AS:i:37 NM:i:11 ZS:i:28
3:153409880:F:224;None;3,153410143,G,A/1    0   chr1    523 4   2M1D32M1D3M1D6M1D8M *   0   0   GAAGAGTTAATTTAAGTCACTTCAAACAGATTACGTATCTTTTTTTTCCCT *   AS:i:42 NM:i:16 ZS:i:41
Y:26750420:R:-132;None;None/1   0   chr1    120 4   2M1I4M3D3M1I7M2I9M2D6M1I8M  *   0   0   AACAACAGAAGTTAATTAGCTTCAAAAATACTTTATATTTGCAA    *   AS:i:32 NM:i:16 ZS:i:29
13:91170622:R:-276;None;None/1  0   chr1    302 4   8M1D8M1D3M2D6M1D4M2I2M1D2M3D5M1I4M  *   0   0   CATTTATTGTTGTTTTTAAAGATTAAATGATTAAATGTTTCAAAA   *   AS:i:32 NM:i:18 ZS:i:30
15:37079528:R:-240;None;None/1  0   chr1    4   5   4M2D4M1D9M1I3M4I16M1I3M1D4M2D5M *   0   0   ACAGTGATGCCAAGCCAGTGGGTTTTAGCTTGTGGAGTTCCATAGGAGCGATGC  *   AS:i:30 NM:i:22 ZS:i:23
9:92308501:R:-176;None;None/1   0   chr1    142 4   4M3I5M4D10M2D4M1I2M2I6M5D1M1D6M2D3M *   0   0   AATAACCATAAAAATGGGCAAAGCAGCCTTCAGGGCTGCTGTTTCTA *   AS:i:26 NM:i:25 ZS:i:26
...
```

What is the output?

Please check the document "The SAM Format Specification" at: http://samtools.github.io/hts-specs/SAMv1.pdf for the full description.

The additional optional field "ZS" indicates the suboptimal alignment score. For example, the 1st record in the upper example means the optimal alignment score of the given sequence is 37; the suboptimal alignment score is 28; the mismatch and INDEL base count within the aligned fragment of the read is 11.

2. An example of the BLAST like output:

```
target_name: chr1
query_name: 6:163296599:F:198;None;None/1
optimal_alignment_score: 37	sub-optimal_alignment_score: 28	strand: +	target_begin: 453	target_end: 492	query_begin: 17	query_end: 51

Target:      453    CCAATGCCACAAAACATCTGTCTCTAACTGGTG--TGTGTGT    492
                    |||  ||| ||||  |||||| | ||| |||||  |*|||||
Query:        17    CCA--GCC-CAAA--ATCTGT-TTTAA-TGGTGGATTTGTGT    51

target_name: chr1
query_name: 3:153409880:F:224;None;3,153410143,G,A/1
optimal_alignment_score: 42	sub-optimal_alignment_score: 41	strand: +	target_begin: 523	target_end: 577	query_begin: 3	query_end: 53

Target:      523    GAGAGAGAAAATTTCACTCCCTCCATAAATCTCACAGTATTCTTTTCTTTTTCCT    577
                    || ||||**|||||*|*||*||*||*|*|**|*|| ||| |||||| ||||*|||
Query:         3    GA-AGAGTTAATTTAAGTCACTTCAAACAGATTAC-GTA-TCTTTT-TTTTCCCT    53
...
```

## Python interface

### How to use the Python wrapper `ssw_lib.py`

`ssw_lib.py` is a Python wrapper that fully supports APIs of the C library. To use this Python library, C programming knowledge is not required.

To use the Python wrapper, please:

1. Compile the `src` folder by either using the `makefile` or by compiling a dynamic shared library with `gcc`:

```
gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h
```

2. Put `libssw.so` and `ssw_lib.py` in the same directory of your own program files or directories in sys.paths.

3. The `LD_LIBRARY_PATH` environment variable may need to be modified to include the directory of the dynamic library `libssw.so` by one of the two following mathods:

    * ```export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_of_libssw.so```
    * For a definitive inclusion edit `/etc/ld.so.conf` and add the path of the `libssw.so`. Then, update the cache by `/sbin/ldconfig`.

4. In a Python script or in a interactive interpreter, import the CSsw class by: `from ssw_lib import CSsw` or `import ssw_lib` and then call `ssw_lib.CSsw`.

5. Call APIs with input parameters and parse the results (Please see `pyssw.py` as an example).

### Run pyssw standalone

```
usage: pyssw.py [-h] [-l SLIBPATH] [-m NMATCH] [-x NMISMATCH] [-o NOPEN]
                [-e NEXT] [-p] [-a SMATRIX] [-c] [-f NTHR] [-r] [-s] [-header]
                [target] [query]

positional arguments:
  target                targe file
  query                 query file

optional arguments:
  -h, --help            show this help message and exit
  -l SLIBPATH, --sLibPath SLIBPATH
                        path of libssw.so
  -m NMATCH, --nMatch NMATCH
                        a positive integer as the score for a match in genome
                        sequence alignment. [default: 2]
  -x NMISMATCH, --nMismatch NMISMATCH
                        a positive integer as the score for a mismatch in
                        genome sequence alignment. [default: 2]
  -o NOPEN, --nOpen NOPEN
                        a positive integer as the penalty for the gap opening
                        in genome sequence alignment. [default: 3]
  -e NEXT, --nExt NEXT  a positive integer as the penalty for the gap
                        extension in genome sequence alignment. [default: 1]
  -p, --bProtien        Do protein sequence alignment. Without this option,
                        the ssw_test will do genome sequence alignment.
                        [default: False]
  -a SMATRIX, --sMatrix SMATRIX
                        a file for either Blosum or Pam weight matrix.
                        [default: Blosum50]
  -c, --bPath           Return the alignment path. [default: False]
  -f NTHR, --nThr NTHR  a positive integer. Only output the alignments with
                        the Smith-Waterman score >= N.
  -r, --bBest           The best alignment will be picked between the original
                        read alignment and the reverse complement read
                        alignment. [default: False]
  -s, --bSam            Output in SAM format. [default: no header]
  -header, --bHeader    If -s is used, include header in SAM output.
```

### pyssw output

The software can output SAM format or BLAST like format results.

### Speed and memory usage of pyssw

The speed and memory are about the same as the C library.

## Java interface

The Java wrapper is a thin JNI (Java Native Interface) wrapper around the native C implementation.

### Building

Only the C, C++, and C shared libraries are generated from the default make goal, and as such, the Java interface must be built explicitly.

1. Ensure `javac` and `jar` are in `PATH`.
2. Ensure `JAVA_HOME` is set to an installed JRE or JDK, or the JNI include directory is included in the C system library search path.
3. `make java` from the src directory.
4. `libsswjni.so` and `ssw.jar` should be built.

### How to use the Java wrapper

The Java wrapper consist of the following components:

* `libsswjni.so`: native C library exposing the JNI entry points
* `ssw.jar` is a Java library containing the Java interface to the native C library. This small wrapper library is composed of:

    * `ssw.Aligner` Java class: a thread-safe static class that exposes two `align()` methods. The first exposes the SSW C library directly. No error checking is performed on arguments passed to this method and misuse is highly likely to crash the JVM. The second `align()` method is a more user-friendly entry point that exposes a simpler API and performs some basic error checking.
    * `ssw.Alignment` Java class: this class stores alignment results. Each each field has a direct correspondence and identical meaning to the C s_align struct.
    * `ssw.Example` Java class: Java version of the example_c sample code. Run `java -jar ssw.jar` to execute the sample.

To use the library, either reference the `ssw.jar` or including the Aligner and Alignment classes directly. As for any JNI library, the native library must be loaded (using `System.loadLibrary("sswjni")` or similar) before invokation of native methods. For the JVM to find the library, ensure that either the library is included in the `LD_LIBRARY_PATH` environment variable, or `-Djava.library.path=<directory containing libsswjni.so>` is set on the Java command line.

## Citation

Please cite this paper, if you need:
https://doi.org/10.1371/journal.pone.0082138
