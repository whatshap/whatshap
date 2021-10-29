import logging
from abc import ABC, abstractmethod
from typing import Dict

from math import log

import networkx as nx

from whatshap.core import Read, ReadSet

logger = logging.getLogger(__name__)


class ReadMergerBase(ABC):
    @abstractmethod
    def merge(self, readset: ReadSet) -> ReadSet:
        pass


class ReadMerger(ReadMergerBase):
    def __init__(
        self, error_rate: float, max_error_rate: float, positive_threshold, negative_threshold
    ):
        """
        error_rate: the probability that a nucleotide is wrong
        max_error_rate: the maximum error rate of any edge of the read
            merging graph allowed before we discard it
        positive_threshold: The threshold of the ratio between the probabilities
            that a pair of reads come from the same haplotype and different
            haplotypes
        negative_threshold: The threshold of the ratio between the probabilities
            that a pair of reads come from the same haplotype and different haplotypes.
        """
        self._error_rate = error_rate
        self._max_error_rate = max_error_rate
        self._positive_threshold = positive_threshold
        self._negative_threshold = negative_threshold

    def merge(self, readset: ReadSet) -> ReadSet:
        """
        Return a set of reads after merging together subsets of reads
        (into super reads) from an input readset according to a
        probabilistic model of how likely sets of reads are to appear
        together on one haplotype and on opposite haplotypes.

        readset -- the input .core.ReadSet object
        """
        logger.info(
            "Merging %d reads with error rate %.2f, maximum error rate %.2f, "
            "positive threshold %d and negative threshold %d ...",
            len(readset),
            self._error_rate,
            self._max_error_rate,
            self._positive_threshold,
            self._negative_threshold,
        )
        logger.debug("Merging started.")
        gblue = nx.Graph()
        gnotblue = nx.Graph()

        # Probability that any nucleotide is wrong
        error_rate = self._error_rate
        logger.debug("Error Rate: %s", error_rate)

        # If an edge has too many errors, we discard it since it is not reliable
        logger.debug("Max Error Rate: %s", self._max_error_rate)

        # Threshold of the ratio between the probabilities that the two reads come from
        # the same side or from different sides
        thr = self._positive_threshold
        logger.debug("Positive Threshold: %s", thr)

        # Threshold_neg is a more conservative threshold for the evidence
        # that two reads should not be clustered together.
        thr_neg = self._negative_threshold
        logger.debug("Negative Threshold: %s", thr_neg)

        thr_diff = 1 + int(log(thr, (1 - error_rate) / (error_rate / 3)))
        thr_neg_diff = 1 + int(log(thr_neg, (1 - error_rate) / (error_rate / 3)))
        logger.debug("Thr. Diff.: %s - Thr. Neg. Diff.: %s", thr_diff, thr_neg_diff)

        logger.debug("Start reading the reads...")
        reads = []
        queue = {}
        for i, read in enumerate(readset):
            alleles = []
            orgn = []
            for variant in read:
                position = variant.position
                allele = variant.allele
                quality = variant.quality

                orgn.append((position, allele, quality))
                assert allele in (0, 1)
                alleles.append(allele)
            reads.append(orgn)

            begin = read[0].position
            end = begin + len(alleles)
            gblue.add_node(i, begin=begin, end=end)
            gnotblue.add_node(i, begin=begin, end=end)
            queue[i] = {"begin": begin, "end": end, "alleles": alleles}
            for x in [id for id in queue.keys() if queue[id]["end"] <= begin]:  # type: ignore
                del queue[x]
            for j in queue.keys():
                if i == j:
                    continue
                match, mismatch = eval_overlap(queue[j], queue[i])
                if (
                    match + mismatch >= thr_neg_diff
                    and min(match, mismatch) / (match + mismatch) <= self._max_error_rate
                    and match - mismatch >= thr_diff
                ):
                    gblue.add_edge(j, i, match=match, mismatch=mismatch)
                    if mismatch - match >= thr_neg_diff:
                        gnotblue.add_edge(j, i, match=match, mismatch=mismatch)

        logger.debug("Finished reading the reads.")
        logger.debug("Number of reads: %s", len(reads))
        logger.debug("Blue Graph")
        logger.debug(
            "Nodes: %s - Edges: %s - ConnComp: %s",
            nx.number_of_nodes(gblue),
            nx.number_of_edges(gblue),
            len(list(nx.connected_components(gblue))),
        )
        logger.debug("Non-Blue Graph")
        logger.debug(
            "Nodes: %s - Edges: %s - ConnComp: %s",
            nx.number_of_nodes(gnotblue),
            nx.number_of_edges(gnotblue),
            len(list(nx.connected_components(gnotblue))),
        )

        # We consider the notblue edges as an evidence that two reads
        # should not be merged together
        # Since we want to merge each blue connected components into
        # a single superread, we check each notblue edge (r1, r2) and
        # we remove some blue edges so that r1 and r2 are not in the
        # same blue connected component

        blue_component = {}
        current_component = 0
        for conncomp in nx.connected_components(gblue):
            for v in conncomp:
                blue_component[v] = current_component
            current_component += 1

        for (u, v) in gnotblue.edges():
            if blue_component[u] != blue_component[v]:
                # Keep only the notblue edges that are inside a blue connected component
                continue
            while v in nx.node_connected_component(gblue, u):
                path = nx.shortest_path(gblue, source=u, target=v)
                # Remove the edge with the smallest support
                # A better strategy is to weight each edge with -log p
                # and remove the minimum (u,v)-cut
                w, x = min(
                    zip(path[:-1], path[1:]),
                    key=lambda p: gblue[p[0]][p[1]]["match"] - gblue[p[0]][p[1]]["mismatch"],
                )
                gblue.remove_edge(w, x)

        # Merge blue components (somehow)
        logger.debug("Started Merging Reads...")
        superreads: Dict = {}  # superreads given by the clusters (if clustering)
        representative = {}  # cluster representative of a read in a cluster

        for cc in nx.connected_components(gblue):
            if len(cc) == 1:
                continue
            r = min(cc)
            superreads[r] = {}
            for i in cc:
                representative[i] = r

        for id in range(len(reads)):
            if id in representative:
                for position, allele, quality in reads[id]:
                    r = representative[id]
                    if position not in superreads[r]:
                        superreads[r][position] = [0, 0]
                    superreads[r][position][allele] += quality

        merged_reads = ReadSet()
        readn = 0
        for id in range(len(reads)):
            read = Read(f"read{readn}")
            readn += 1
            if id in representative:
                if id == representative[id]:
                    for position in sorted(superreads[id]):
                        z = superreads[id][position]
                        allele = 0 if z[0] >= z[1] else 1
                        read.add_variant(position, allele, abs(z[1] - z[0]))
                    merged_reads.add(read)
            else:
                for position, allele, quality in reads[id]:
                    read.add_variant(position, allele, quality)
                merged_reads.add(read)

        logger.debug("Finished merging reads.")
        logger.info(
            "... after merging: merged %d reads into %d reads", len(readset), len(merged_reads)
        )

        return merged_reads


class DoNothingReadMerger(ReadMergerBase):
    def merge(self, readset):
        return readset


def eval_overlap(n1, n2):
    """
    Return a tuple containing the number of matches (resp.,
    mismatches) between a pair (n1,n2) of overlapping reads
    """
    hang1 = n2["begin"] - n1["begin"]
    overlap = zip(n1["alleles"][hang1:], n2["alleles"])
    match = mismatch = 0
    for (c1, c2) in overlap:
        if c1 == c2:
            match += 1
        else:
            mismatch += 1
    return match, mismatch
