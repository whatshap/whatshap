#!/usr/bin/env python

from argparse import ArgumentParser
import sys
import os
import subprocess


__author__ = "Amal Zine el Aabidine, Jasmijn Baaijens and Xiao Luo"

usage = """

Convert SFO output to overlaps format for SAVAGE. This script also uses the read
pairing-information to output paired-end read overlaps if present.

"""

def main():
    parser = ArgumentParser(description=usage)
    parser.add_argument('--in', dest='infile', type=str, required=True, help="input file: overlap results from SFO")
    parser.add_argument('--out', dest='outfile', type=str, required=True, help="output file: overlap results for SAVAGE")
    parser.add_argument('--num_singles', dest='num_singles', type=int, required=True, help="number of single-end reads")
    parser.add_argument('--num_pairs', dest='num_pairs', type=int, required=True, help="number of paired-end reads")
    args = parser.parse_args()

    # read SFO results and add two fields for the original read IDs
    # (the two ends of a paired-end read have the same original ID)
    # write the result to a temporary overlaps file:
    outdir = os.path.dirname(os.path.abspath(args.infile))
    tmp_file = outdir+'/tmp_overlaps.txt'
    with open(args.infile, 'r') as f_in:
        with open(tmp_file, 'w') as f_tmp:
            for line in f_in:
                sfo_line = line.strip('\n').split()
                assert len(sfo_line) == 8
                # [sfo_idA, sfo_idB, ori, OHA, OHB, OLA, OLB, K] = sfo_line
                idA = int(sfo_line[0])
                idB = int(sfo_line[1])
                new_idA = get_original_id(idA, args.num_singles, args.num_pairs)
                new_idB = get_original_id(idB, args.num_singles, args.num_pairs)
                if new_idA > new_idB:
                    if sfo_line[2] == 'I':
                        flipped_tup = flip_I(sfo_line)
                    else:
                        flipped_tup = flip_N(sfo_line)
                    flipped_line = '\t'.join(flipped_tup) + '\n'
                    newline = str(new_idB) + "\t" + str(new_idA) + "\t" + flipped_line
                else:
                    newline = str(new_idA) + "\t" + str(new_idB) + "\t" + line
                f_tmp.write(newline)

    # sort the resulting overlaps file
    subprocess.check_call('sort -k1,1n -k2,2n -k3,3n -k4,4n %s | uniq > %s/sorted_overlaps.txt' % (tmp_file,outdir), shell=True)
    subprocess.check_call('mv %s/sorted_overlaps.txt %s' % (outdir,tmp_file), shell=True)

    # read the sorted overlaps and process line by line
    s_s_count = 0
    p_count = 0
    tmp2_file = outdir+ '/tmp_overlaps_paired.txt'
    with open(tmp_file, 'r') as f_tmp:
        with open(tmp2_file, 'w') as f_out:
            candidates = []
            for line in f_tmp:
                sfo_line = line.strip('\n').split()
                assert len(sfo_line) == 10
                # [idA, idB, sfo_idA, sfo_idB, ori, OHA, OHB, OLA, OLB, K] = sfo_line
                idA = int(sfo_line[0])
                idB = int(sfo_line[1])
                if idA == idB: # self-overlap
                    continue

                OLA = int(sfo_line[7]) # read A bases inside overlap
                OLB = int(sfo_line[8]) # read B bases inside overlap
#                if OLA != OLB: # indel in overlap
#                    continue

                is_paired_A = is_paired(idA, args.num_singles, args.num_pairs)
                is_paired_B = is_paired(idB, args.num_singles, args.num_pairs)
                if (not is_paired_A) and (not is_paired_B):
                    # single-single overlap
                    s_s_overlap = get_s_s_overlap(sfo_line)
                    if len(s_s_overlap) > 0:
                        overlap_line = '\t'.join(s_s_overlap) + '\n'
                        f_out.write(overlap_line)
                        s_s_count += 1
                else:
                    # paired-end read involved in overlap
                    current_id_pair = [str(idA), str(idB)]
                    if len(candidates) > 0:
                        candidates_ids = candidates[0][0:2]
                        assert idA >= int(candidates_ids[0])
                        if idA == int(candidates_ids[0]):
                            assert idB >= int(candidates_ids[1])
                        if candidates_ids != current_id_pair:
                            paired_overlaps = match_candidates(candidates, is_paired_A, is_paired_B)
                            for p_overlap in paired_overlaps:
                                if len(p_overlap) > 0:
                                    overlap_line = '\t'.join(p_overlap) + '\n'
                                    f_out.write(overlap_line)
                                    p_count += 1
                            candidates = []
                    candidates.append(sfo_line)
    # print("total overlap count: %d" % (s_s_count + p_count))
    # print("of which single-single: %d" % s_s_count)

    # remove duplicate overlaps from the resulting overlaps file
    subprocess.check_call("uniq %s > %s" % (tmp2_file, args.outfile), shell=True)
    # remove temporary overlaps file
    subprocess.check_call("rm %s %s" % (tmp_file, tmp2_file), shell=True)


def flip_N(sfo_tup):
    # sfo_tuple = [sfo_idA, sfo_idB, ori, OHA, OHB, OLA, OLB, K]
    OHA = -1*int(sfo_tup[3])
    OHB = -1*int(sfo_tup[4])
    flipped = [sfo_tup[1], sfo_tup[0], sfo_tup[2], str(OHA), str(OHB), sfo_tup[6], sfo_tup[5], sfo_tup[7]]
    return flipped

def flip_I(sfo_tup):
    # sfo_tuple = [sfo_idA, sfo_idB, ori, OHA, OHB, OLA, OLB, K]
    flipped = [sfo_tup[1], sfo_tup[0], sfo_tup[2], sfo_tup[4], sfo_tup[3], sfo_tup[6], sfo_tup[5], sfo_tup[7]]
    return flipped

def is_paired(ID, num_singles, num_pairs):
    # check if a read is paired based on its read ID
    if num_pairs == 0:
        paired = False
    else:
        assert ID >= 0 and ID < num_singles + num_pairs
        if ID >= num_singles:
            paired = True
        else:
            paired = False
    return paired

def get_original_id(sfo_ID, num_singles, num_pairs):
    # translate the read ID used by SFO to the original read ID
    if num_pairs == 0:
        return sfo_ID
    assert sfo_ID >= 0 and sfo_ID < num_singles + 2*num_pairs
    if sfo_ID < num_singles + num_pairs:
        # single-end read or /1 read of a pair
        original_ID = sfo_ID
    else:
        # /2 read of a pair
        original_ID = sfo_ID - num_pairs
    return original_ID


def get_s_s_overlap(sfo_line):
    # translate SFO style overlap to SAVAGE style overlap
    # sfo format: [idA, idB, sfo_idA, sfo_idB, ori, OHA, OHB, OLA, OLB, K]
    # savage format: ID1 ID2 POS1 POS2 ORD ORI1 ORI2 PERC1 PERC2 LEN1 LEN2 TYPE1 TYPE2
    idA = sfo_line[0]
    idB = sfo_line[1]
    OHA = int(sfo_line[5]) # read A bases outside overlap
    OHB = int(sfo_line[6]) # read B bases outside overlap
    OLA = int(sfo_line[7]) # read A bases inside overlap
    OLB = int(sfo_line[8]) # read B bases inside overlap
    ori = "+" if sfo_line[4] == "N" else "-"
    ovlen = min(OLA, OLB)
    if OHA >= 0: # read A is first
        if OHB >= 0:
            readlenA = OLA + OHA
            readlenB = OLB + OHB
        else:
            readlenA = OLA + OHA + -OHB
            readlenB = OLB
        id1 = idA
        id2 = idB
        pos1 = str(OHA)
        pos2 = "-"
        order = "-"
        ori1 = "+"
        ori2 = ori
    else: # read B is first
        if OHB >= 0:
            readlenA = OLA
            readlenB = -OHA + OLB + OHB
        else:
            readlenA = OLA + -OHB
            readlenB = -OHA + OLB
        id1 = idB
        id2 = idA
        pos1 = str(-1*OHA)
        pos2 = "-"
        order = "-"
        ori1 = ori
        ori2 = "+"
    minreadlen = min(readlenA, readlenB)
    perc = min(round(100*ovlen/minreadlen), 100)
    perc1 = "{:.0f}".format(perc)
    perc2 = "-"
    len1 = str(ovlen)
    len2 = "-"
    type1 = "s"
    type2 = "s"
    assert minreadlen > 0
    overlap = [id1, id2, pos1, pos2, order, ori1, ori2, perc1, perc2, len1, len2, type1, type2]
    return overlap


def match_candidates(candidates, typeA, typeB):
    # using maxoverlaps for SFO should result in only one candidate per
    # read end, so at most two candidates for every paired-end read;
    # for a paired end overlap, we need exactly two candidates.
    overlaps = []
    if len(candidates) < 2:
        return []
    if len(candidates) >= 2:
#        print candidates
        for i in xrange(len(candidates)):
            cand1 = candidates[i]
            for j in range(i+1, len(candidates)):
                cand2 = candidates[j]
                overlap = find_paired_overlap(cand1, cand2, typeA, typeB)
                if len(overlap) > 0:
                    overlaps.append(overlap)
    return overlaps


def find_paired_overlap(cand1, cand2, typeA, typeB):
#    assert len(candidates) == 2
#    cand1 = candidates[0]
#    cand2 = candidates[1]
    overlap1 = []
    overlap2 = []
    # sfo format: [idA, idB, sfo_idA, sfo_idB, ori, OHA, OHB, OLA, OLB, K]
    if cand1[4] != cand2[4]:
        # orientations don't match; since SFO outputs the orientation of
        # read B, which is always the read with largest ID, the orientations
        # of both read ends should agree
        return []
    cand1_id1 = int(cand1[2])
    cand1_id2 = int(cand1[3])
    cand2_id1 = int(cand2[2])
    cand2_id2 = int(cand2[3])
    if typeA and typeB:
        # paired-paired
        if cand1[4] == 'N':
            # forward-forward
            if cand1_id1 < cand2_id1 and cand1_id2 < cand2_id2:
                overlap1 = get_s_s_overlap(cand1)
                overlap2 = get_s_s_overlap(cand2)
            elif cand1_id1 > cand2_id1 and cand1_id2 > cand2_id2:
                overlap1 = get_s_s_overlap(cand2)
                overlap2 = get_s_s_overlap(cand1)
        elif cand1[4] == 'I':
            # forward-reverse
            if cand1_id1 < cand2_id1 and cand1_id2 > cand2_id2:
                overlap1 = get_s_s_overlap(cand1)
                overlap2 = get_s_s_overlap(cand2)
            elif cand1_id1 > cand2_id1 and cand1_id2 < cand2_id2:
                overlap1 = get_s_s_overlap(cand2)
                overlap2 = get_s_s_overlap(cand1)
    elif typeA and not typeB:
        # paired-single
        cand1_pos1 = int(cand1[5])
        cand2_pos1 = int(cand2[5])
        if cand1[4] == 'N':
            # forward-forward
            if (cand1_id1 < cand2_id1 and cand1_pos1 < cand2_pos1):
                overlap1 = get_s_s_overlap(cand1)
                overlap2 = get_s_s_overlap(cand2)
            elif (cand1_id1 > cand2_id1 and cand1_pos1 > cand2_pos1):
                overlap1 = get_s_s_overlap(cand2)
                overlap2 = get_s_s_overlap(cand1)
        elif cand1[4] == 'I':
            # forward-reverse
            if (cand1_id1 < cand2_id1 and cand1_pos1 > cand2_pos1):
                overlap1 = get_s_s_overlap(cand2)
                overlap2 = get_s_s_overlap(cand1)
            elif (cand1_id1 > cand2_id1 and cand1_pos1 < cand2_pos1):
                overlap1 = get_s_s_overlap(cand1)
                overlap2 = get_s_s_overlap(cand2)
    else:
        # single-paired
        cand1_pos1 = int(cand1[5])
        cand2_pos1 = int(cand2[5])
        if cand1[4] == 'N':
            # forward-forward
            if (cand1_id2 < cand2_id2 and cand1_pos1 < cand2_pos1):
                overlap1 = get_s_s_overlap(cand1)
                overlap2 = get_s_s_overlap(cand2)
            elif (cand1_id2 > cand2_id2 and cand1_pos1 > cand2_pos1):
                overlap1 = get_s_s_overlap(cand2)
                overlap2 = get_s_s_overlap(cand1)
        elif cand1[4] == 'I':
            # forward-reverse
            if (cand1_id2 < cand2_id2 and cand1_pos1 > cand2_pos1):
                overlap1 = get_s_s_overlap(cand2)
                overlap2 = get_s_s_overlap(cand1)
            elif (cand1_id2 > cand2_id2 and cand1_pos1 < cand2_pos1):
                overlap1 = get_s_s_overlap(cand1)
                overlap2 = get_s_s_overlap(cand2)

    if len(overlap1) > 0 and len(overlap2) > 0:
        # now combine both overlaps into one paired-end overlap
        if overlap1[0] == cand1[0]:
            assert overlap1[1] == cand1[1]
            type1 = "p" if typeA else "s"
            type2 = "p" if typeB else "s"
        else:
            assert overlap1[1] == cand1[0]
            assert overlap1[0] == cand1[1]
            type1 = "p" if typeB else "s"
            type2 = "p" if typeA else "s"
        final_overlap = merge_overlaps(overlap1, overlap2, type1, type2)
        return final_overlap
    else:
        return []

def merge_overlaps(overlap1, overlap2, type1, type2):
    # merge two matching single-single overlaps into a paired-end overlap
    # savage format: ID1 ID2 POS1 POS2 ORD ORI1 ORI2 PERC1 PERC2 LEN1 LEN2 TYPE1 TYPE2
    overlap = overlap1
    overlap[11] = type1
    overlap[12] = type2
    if type1 == "p" and type2 == "p":
        # take care of ord
        if overlap1[0] != overlap2[0]:
            assert overlap1[0] == overlap2[1]
            overlap[4] = "2"
        else:
            overlap[4] = "1"
    overlap[3] = overlap2[2] # pos2
    overlap[8] = overlap2[7] # perc2
    overlap[10] = overlap2[9] # len2
    return overlap


if __name__ == '__main__':
    sys.exit(main())
