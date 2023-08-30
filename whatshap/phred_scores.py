"""
Takes as input the chromosome wise ref-read kmer pair counts from whatshap learn,
calculates the probability for observing each ref-read kmer pair across the whole genome
and outputs a table with corresponding phred-scores
"""

import argparse
import math
from pathlib import Path
import csv


def phred(inputfile, outputfile, epsilon_value, kmer_size):
    probs = {}
    probabilities = {}
    seen_comb = {}
    comb_count = {}
    k = int(kmer_size)
    total_comb = float(4**k)
    epsilon = float(epsilon_value)
    path = Path(inputfile)
    glob_path = path.glob("*.txt")
    writer = open(outputfile, "w")
    # Iterate over each chromosome file
    for file in glob_path:
        counts_file = open(file, "r")
        reader = csv.reader(counts_file, delimiter="\t")
        for line in reader:
            ref = line[1]
            read = line[2]
            count = int(line[3])
            new_key = (ref, read)
            if new_key in probs:
                probs[new_key] += count
            else:
                probs[new_key] = count

    for every_key in probs:
        if every_key[0] in seen_comb:
            seen_comb[every_key[0]] += 1
        else:
            seen_comb[every_key[0]] = 1
        if every_key[0] in comb_count:
            comb_count[every_key[0]] += probs[every_key]
        else:
            comb_count[every_key[0]] = probs[every_key]

    for s in probs:
        probability = probs[s] / (comb_count[s[0]] + ((total_comb - seen_comb[s[0]]) * epsilon))
        e_probability = epsilon / (comb_count[s[0]] + ((total_comb - seen_comb[s[0]]) * epsilon))
        phred_score = -10 * math.log10(float(probability))
        e_phred_score = -10 * math.log10(float(e_probability))
        if s[0] in probabilities:
            probabilities[s[0]].append("1")
            writer.write(str(s[0]) + "\t" + str(s[1]) + "\t" + str(phred_score) + "\n")
        else:
            probabilities[s[0]] = list()
            probabilities[s[0]].append("1")
            writer.write(str(s[0]) + "\t" + str(s[1]) + "\t" + str(phred_score) + "\n")
            writer.write(str(s[0]) + "\t" + "-5" + "\t" + str(e_phred_score) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--inputfile",
        type=str,
        help="The path to folder containing ref-read kmer pair counts obtained from whatshap call",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--outputfile",
        type=str,
        help="A table of phred scores based on probilities calculated using kmer pair counts",
        required=True,
    )
    parser.add_argument(
        "-e",
        "--epsilon_value",
        type=str,
        help="The pseudocount value for kmer pairs that are not observed ",
        required=True,
    )
    parser.add_argument("-k", "--kmer_size", type=str, help="kmer length", required=True)
    args = parser.parse_args()

    phred(args.inputfile, args.outputfile, args.epsilon_value, args.kmer_size)
