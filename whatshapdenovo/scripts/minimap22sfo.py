#!/usr/bin/env python
from argparse import ArgumentParser
import sys


__author__ = "Jasmijn Baaijens"
__license__ = "GPL"

usage = """

Convert Pairwise mApping Format (PAF) overlaps to SFO format

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--in', dest='infile', type=str)
    parser.add_argument('--out', dest='outfile', type=str)
    parser.add_argument('-m', dest='min_overlap_len', default=0, type=int)
    parser.add_argument('-p', dest='min_pident', default=0, type=float)
    args = parser.parse_args()

    if not (args.infile and args.outfile):
        print("Specify input and output files.")
        parser.print_help()
        sys.exit(1)

    with open(args.outfile, 'w') as f1:
        with open(args.infile, 'r') as f2:
            too_short_count = 0
            too_div_count = 0
            overlap_count = 0
            for line in f2:
                [qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length] = line.strip('\n').split('\t')[:11]
                #                [qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, qlen, slen] = line.strip('\n').split('\t')
                #                qori = int(qstart) <= int(qend)
                #                sori = int(sstart) <= int(send)
                #                assert qori # otherwise script needs extra cases to find correct pos1
                if int(length) < args.min_overlap_len:
                    too_short_count += 1
                    continue
                if int(matchcount)/float(length) < args.min_pident/100.0:
                    #                    print matchcount, length
                    too_div_count += 1
                    continue

                idA = qseqid
                idB = sseqid
                ori = 'N' if qori == '+' else 'I'
                if ori == 'N':
                    OHA = int(qstart) - int(sstart)
                    OHB = int(slen) - int(sstart) - (int(qlen) - int(qstart))
                else:
                    OHA = int(qstart) - (int(slen) - int(send))
                    OHB = int(send) - (int(qlen) - int(qstart))
                if OHA >= 0:
                    OLA = min(int(qlen) - OHA, int(slen))
                else:
                    OLA = min(int(slen) + OHA, int(qlen))
                OLB = OLA
                #                if int(idA) > int(idB):
                if idA > idB:
                    # swap order such that id1 < id2
                    idA, idB = idB, idA
                    if ori == 'N':
                        OHA *= -1
                        OHB *= -1
                    else:
                        # swap orientations such that id1 sequence is forward
                        OHA, OHB = OHB, OHA
                mismatch = int(length) - int(matchcount)
                assert mismatch >= 0
                sfo_line = '\t'.join([idA, idB, ori, str(OHA), str(OHB), str(OLA), str(OLB), str(mismatch)]) + '\n'
                f1.write(sfo_line)
                overlap_count += 1
            # print("overlaps shorter than %d bp: %d" %(args.min_overlap_len, too_short_count))
            # print("overlaps with less than %d percent identity: %d" %(args.min_pident, too_div_count))
            # print("total overlaps found: ", overlap_count)

if __name__ == '__main__':
    sys.exit(main())
