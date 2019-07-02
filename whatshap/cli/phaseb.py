"""
Phase reads mapped to bubble chains
"""

import logging
import os
import sys
from .core import ReadSet, Read, Pedigree, PedigreeDPTable, NumericSampleIds
from .pedigree import uniform_recombination_map
from .vg_pb2 import Alignment
import stream

__author__ = "Fawaz Dabbaghie"

logger = logging.getLogger(__name__)


class Node:
    def __init__(self, identifier):
        self.id = identifier
        self.seq = ""
        self.seq_len = 0
        self.start = []
        self.end = []
        self.visited = False
        self.which_chain = 0
        self.which_sb = 0
        self.which_b = 0
        self.which_allele = -1


def add_arguments(parser):
    add = parser.add_argument
    add('gfa_file',metavar='GFA', type=str, help='Give the modified GFA file path')
    add('gam_file',metavar='GAM', type=str, help='Give the alignment GAM file path')


def validate(args, parser):
    if not os.path.exists(args.gfa_file):
        parser.error("The GFA file was not found")

    if not os.path.exists(args.gam_file):
        parser.error("The GAM file was not found")


def read_gfa(gfa_file_path, modified=False):
    """
    :param gfa_file_path: gfa graph file.
    :param modified: if I'm reading my modified GFA with extra information for the nodes
    :return: Dictionary of node ids and Node objects.
    """
    if not os.path.exists(gfa_file_path):
        print("the gfa file path you gave does not exists, please try again!")
        sys.exit()

    nodes = dict()
    edges = []
    with open(gfa_file_path, "r") as lines:
        for line in lines:
            if line.startswith("S"):
                if modified:
                    line = line.split("\t")
                    n_id = int(line[1])
                    nodes[n_id] = Node(n_id)
                    nodes[n_id].seq_len = len(line[2])
                    nodes[n_id].seq = str(line[2])
                    # this is my extra columns
                    # constructed as which_chain:which_sb:which_b:which_allele
                    # e.g. 5:0:3:1 (chain 5, bubble 3, allele 1)
                    specifications = str(line[3])
                    specifications = specifications.split(":")
                    nodes[n_id].which_chain = int(specifications[0])
                    nodes[n_id].which_sb = int(specifications[1])
                    nodes[n_id].which_b = int(specifications[2])
                    nodes[n_id].which_allele = int(specifications[3])

                else:

                    line = line.split()
                    n_id = int(line[1])
                    n_len = len(line[2])
                    nodes[n_id] = Node(n_id)
                    nodes[n_id].seq_len = n_len
                    nodes[n_id].seq = str(line[2])

            elif line.startswith("L"):
                edges.append(line)

    for e in edges:
        line = e.split()

        k = int(line[1])
        neighbor = int(line[3])
        if line[2] == "-":
            from_start = True
        else:
            from_start = False

        if line[4] == "-":
            to_end = True
        else:
            to_end = False

        if from_start is True and to_end is True:  # from start to end L x - y -
            if (neighbor, 1) not in nodes[k].start:
                nodes[k].start.append((neighbor, 1))
            if (k, 0) not in nodes[neighbor].end:
                nodes[neighbor].end.append((k, 0))

        elif from_start is True and to_end is False:  # from start to start L x - y +

            if (neighbor, 0) not in nodes[k].start:
                nodes[k].start.append((neighbor, 0))

            if (k, 0) not in nodes[neighbor].start:
                nodes[neighbor].start.append((k, 0))

        elif from_start is False and to_end is False:  # from end to start L x + y +
            if (neighbor, 0) not in nodes[k].end:
                nodes[k].end.append((neighbor, 0))

            if (k, 1) not in nodes[neighbor].start:
                nodes[neighbor].start.append((k, 1))

        elif from_start is False and to_end is True:  # from end to end L x + y -
            if (neighbor, 1) not in nodes[k].end:
                nodes[k].end.append((neighbor, 1))

            if (k, 1) not in nodes[neighbor].end:
                nodes[neighbor].end.append((k, 1))

    return nodes


def build_readsets(nodes, gam_file_path):
    read_sets = {}
    with stream.open(str(gam_file_path), "rb") as instream:
        for data in instream:
            g = Alignment()
            g.ParseFromString(data)

            # construct read object
            read = Read(g.name)  # leave other values as default

            for m in g.path.mapping:
                n_id = int(m.position.node_id)

                if n_id not in nodes:
                    continue

                if nodes[n_id].which_chain not in read_sets:
                    read_sets[nodes[n_id].which_chain] = ReadSet()

                if 0 >= nodes[n_id].which_allel <= 1:  # otherwise it's an end node or a superbubble node
                    # bubbles are already sorted
                    read.add_variant(nodes[n_id].which_b, nodes[n_id].which_allele, 30)

            read_sets[nodes[n_id].which_chain].add_read(read)

    return read_sets


def run_phaseb(gfa_file, gam_file):

    nodes = read_gfa(gfa_file, modified=True)

    all_read_sets = build_readsets(nodes, gam_file)

    for readset in all_read_sets:
        # sort by starting positions of the reads
        readset.sort()

        # positions and recombination costs
        positions = readset.get_positions()
        recombination_cost = uniform_recombination_map(1.26, positions)

        # construct pedigree
        pedigree = Pedigree(NumericSampleIds())

        # assume all genotypes are heterozygous
        distrust_genotypes = False
        genotypes = [1] * len(positions)
        genotype_likelihoods = None
        sample_name = gfa_file.split(".")[0]
        pedigree.add_individual(sample_name, genotypes, genotype_likelihoods)

        # solve MEC
        dp_table = PedigreeDPTable(readset, 1.26, pedigree, distrust_genotypes, positions)

        # cost = dp_table.get_optimal_cost()
        # vector indicating for each read in which bartition (0 or 1) it ended up (ordered according to ReadSet)
        read_partitioning = dp_table.get_optimal_partitioning()
        # overall_components = find_components(positions, readset)

        ## how to get the order of reads in readset:
        ordered_readnames = [read.name for read in readset]

        assert len(read_partitioning) == len(ordered_readnames)


def main(args):
    # here all the main stuff going to happen
    run_phaseb(**vars(args))
