#!/usr/bin/env python
from __future__ import division
from argparse import ArgumentParser
import os
import sys
import random
import gzip

__author__ = "Vincent Luo"

usage = """%prog [options]

Choose a set number of reads randomly to get a sub fastq/fasta file.
"""


def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument("--input", dest="input", type=str, required=True, help="input fastq/fasta filename")
    parser.add_argument("--output", dest="output", type=str, required=True, help="prefix for output files")
    parser.add_argument("--mode", dest="mode", type=str, required=True, help="fastq or fasta")
    parser.add_argument("--num", dest="num", type=int, required=True, help="number of reads/pairs in output files")
    parser.add_argument("--pair", dest="pair", type=bool, required=False, help="if pair end reads, default: False")
    args = parser.parse_args()
    infile = args.input
    mode=args.mode

    if infile.endswith('.gz'):
        n_reads = int(os.popen("zcat %s|wc -l" % (infile)).read().split()[0])
        # n_reads = int(os.popen("zless %s|wc -l" % (fastq)).read().split()[0]) #zcat does not work on Mac
    else:
        n_reads = int(os.popen("wc -l %s" % (infile)).read().split()[0])

    k=0
    if mode == 'fastq':
        k = 8 if (args.pair) else 4
    elif mode=='fasta':
        k = 4 if (args.pair) else 2
    else:
        print('ERROR: mode must be fastq or fasta')
        sys.exit(1)

    n_reads = int(n_reads / k)

    id = [i for i in range(n_reads)]
    id = {((i - 1) * k + 1): 1 for i in random.sample(id, args.num)}

    fw = open("%s.sub.fastq" % args.output, 'w')

    if infile.endswith('.gz'):
        fq=gzip.open(infile, 'rb')
    else:
        fq=open(infile, 'r')

    i = 0
    flag = 0
    for line in fq:
        i += 1
        if i in id.keys() or (flag > 0 and flag < k):
            if infile.endswith('.gz'):
                fw.write(line.decode().strip()+'\n')
            else:
                fw.write(line.strip()+'\n')
            flag += 1
        if flag == k: flag = 0

    fq.close()
    fw.close()


if __name__ == '__main__':
    sys.exit(main())
