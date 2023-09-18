"""
Takes as input the chromosome wise ref-read kmer pair counts from whatshap learn,
calculates the probability for observing each ref-read kmer pair across the whole genome
and outputs a table with corresponding phred-scores
"""

import argparse
import math
from pathlib import Path
import csv
from collections import defaultdict


def phred(inputfile, outputfile, epsilon_value, kmer_size):
    counts = defaultdict(int)
    seen_comb = defaultdict(int)
    comb_count = defaultdict(int)
    probabilities = {}
    k = int(kmer_size)
    total_comb = float(4**k)
    epsilon = float(epsilon_value)
    path = Path(inputfile)
    glob_path = path.glob("*.txt")
    # Iterate over each chromosome file
    for file in glob_path:
        with open(file, "r") as counts_file:
            reader = csv.reader(counts_file, delimiter="\t")
            for line in reader:
                ref = line[1]
                read = line[2]
                count = int(line[3])
                counts[(ref, read)] += count

    for every_key in counts:
        seen_comb[every_key[0]] += 1
        comb_count[every_key[0]] += counts[every_key]

    with open(outputfile, "w") as writer:
        for s in counts:
            probability = counts[s] / (
                comb_count[s[0]] + ((total_comb - seen_comb[s[0]]) * epsilon)
            )
            e_probability = epsilon / (
                comb_count[s[0]] + ((total_comb - seen_comb[s[0]]) * epsilon)
            )
            phred_score = -10 * math.log10(float(probability))
            e_phred_score = -10 * math.log10(float(e_probability))
            if s[0] not in probabilities:
                probabilities[s[0]] = 1
                print(s[0], -5, e_phred_score, sep="\t", file=writer)
            print(s[0], s[1], phred_score, sep="\t", file=writer)


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
