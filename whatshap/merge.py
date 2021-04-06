import logging
from math import log

from networkx import (
    Graph,
    number_of_nodes,
    number_of_edges,
    connected_components,
    node_connected_component,
    shortest_path,
)
from whatshap.core import Read, ReadSet

logger = logging.getLogger(__name__)


class ReadMerger:
    def __init__(self, error_rate, max_error_rate, positive_threshold, negative_threshold):
        self._error_rate = error_rate
        self._max_error_rate = max_error_rate
        self._positive_threshold = positive_threshold
        self._negative_threshold = negative_threshold

    def merge(self, readset):
        """
        Return a set of reads after merging together subsets of reads
        (into super reads) from an input readset according to a
        probabilistic model of how likely sets of reads are to appear
        together on one haplotype and on opposite haplotypes.
        readset -- the input .core.ReadSet object
        error_rate -- the probability that a nucleotide is wrong
        max_error_rate -- the maximum error rate of any edge of the read
        merging graph allowed before we discard it
        threshold -- the threshold of the ratio between the probabilities
        that a pair ' 'of reads come from the same haplotype and different
        haplotypes
        neg_threshold -- The threshold of the ratio between the
        probabilities that a pair of reads come from the same haplotype
        and different haplotypes.
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
        gblue = Graph()
        gred = Graph()
        gnotblue = Graph()
        gnotred = Graph()

        # Probability that any nucleotide is wrong
        error_rate = self._error_rate
        logger.debug("Error Rate: %s", error_rate)

        # If an edge has too many errors, we discard it since it is not reliable
        max_error_rate = self._max_error_rate
        logger.debug("Max Error Rate: %s", max_error_rate)

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
        id = 0
        orig_reads = {}
        queue = {}
        reads = {}
        for read in readset:
            id += 1
            begin_str = read[0][0]
            snps = []
            orgn = []
            for variant in read:

                site = variant[0]
                zyg = variant[1]
                qual = variant[2]

                orgn.append([str(site), str(zyg), str(qual)])
                if int(zyg) == 0:
                    snps.append("G")
                else:
                    snps.append("C")

            begin = int(begin_str)
            end = begin + len(snps)
            orig_reads[id] = orgn

            gblue.add_node(id, begin=begin, end=end, sites="".join(snps))
            gnotblue.add_node(id, begin=begin, end=end, sites="".join(snps))
            gred.add_node(id, begin=begin, end=end, sites="".join(snps))
            gnotred.add_node(id, begin=begin, end=end, sites="".join(snps))
            queue[id] = {"begin": begin, "end": end, "sites": snps}
            reads[id] = {"begin": begin, "end": end, "sites": snps}
            for x in [id for id in queue.keys() if queue[id]["end"] <= begin]:
                del queue[x]
            for id1 in queue.keys():
                if id == id1:
                    continue
                match, mismatch = eval_overlap(queue[id1], queue[id])
                if (
                    match + mismatch >= thr_neg_diff
                    and min(match, mismatch) / (match + mismatch) <= max_error_rate
                    and match - mismatch >= thr_diff
                ):
                    gblue.add_edge(id1, id, match=match, mismatch=mismatch)
                    if mismatch - match >= thr_diff:
                        gred.add_edge(id1, id, match=match, mismatch=mismatch)
                    if match - mismatch >= thr_neg_diff:
                        gnotred.add_edge(id1, id, match=match, mismatch=mismatch)
                    if mismatch - match >= thr_neg_diff:
                        gnotblue.add_edge(id1, id, match=match, mismatch=mismatch)

        logger.debug("Finished reading the reads.")
        logger.debug("Number of reads: %s", id)
        logger.debug("Blue Graph")
        logger.debug(
            "Nodes: %s - Edges: %s - ConnComp: %s",
            number_of_nodes(gblue),
            number_of_edges(gblue),
            len(list(connected_components(gblue))),
        )
        logger.debug("Non-Blue Graph")
        logger.debug(
            "Nodes: %s - Edges: %s - ConnComp: %s",
            number_of_nodes(gnotblue),
            number_of_edges(gnotblue),
            len(list(connected_components(gnotblue))),
        )
        logger.debug("Red Graph")
        logger.debug(
            "Nodes: %s - Edges: %s - ConnComp: %s",
            number_of_nodes(gred),
            number_of_edges(gred),
            len(list(connected_components(gred))),
        )
        logger.debug("Non-Red Graph")
        logger.debug(
            "Nodes: %s - Edges: %s - ConnComp: %s",
            number_of_nodes(gnotred),
            number_of_edges(gnotred),
            len(list(connected_components(gnotred))),
        )

        # We consider the notblue edges as an evidence that two reads
        # should not be merged together
        # Since we want to merge each blue connected components into
        # a single superread, we check each notblue edge (r1, r2) and
        # we remove some blue edges so that r1 and r2 are not in the
        # same blue connected component

        blue_component = {}
        current_component = 0
        for conncomp in connected_components(gblue):
            for v in conncomp:
                blue_component[v] = current_component
            current_component += 1

        # Keep only the notblue edges that are inside a blue connected component
        good_notblue_edges = [
            (v, w) for (v, w) in gnotblue.edges() if blue_component[v] == blue_component[w]
        ]

        for (u, v) in good_notblue_edges:
            while v in node_connected_component(gblue, u):
                path = shortest_path(gblue, source=u, target=v)
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
        superreads = {}  # superreads given by the clusters (if clustering)
        rep = {}  # cluster representative of a read in a cluster

        for cc in connected_components(gblue):
            if len(cc) > 1:
                r = min(cc)
                superreads[r] = {}
                for id in cc:
                    rep[id] = r

        for id in orig_reads:
            if id in rep:
                for tok in orig_reads[id]:
                    site = int(tok[0])
                    zyg = int(tok[1])
                    qual = int(tok[2])
                    r = rep[id]
                    if site not in superreads[r]:
                        superreads[r][site] = [0, 0]
                    superreads[r][site][zyg] += qual

            merged_reads = ReadSet()
            readn = 0
            for id in orig_reads:
                read = Read("read" + str(readn))
                readn += 1
                if id in rep:
                    if id == rep[id]:
                        for site in sorted(superreads[id]):
                            z = superreads[id][site]
                            if z[0] >= z[1]:
                                read.add_variant(site, 0, z[0] - z[1])

                            elif z[1] > z[0]:
                                read.add_variant(site, 1, z[1] - z[0])
                        merged_reads.add(read)
                else:
                    for tok in orig_reads[id]:
                        read.add_variant(int(tok[0]), int(tok[1]), int(tok[2]))
                    merged_reads.add(read)

        logger.debug("Finished merging reads.")
        logger.info(
            "... after merging: merged %d reads into %d reads", len(readset), len(merged_reads)
        )

        return merged_reads


class DoNothingReadMerger:
    @staticmethod
    def merge(readset):
        return readset


def eval_overlap(n1, n2):
    """
    Return a tuple containing the number of matches (resp.,
    mismatches) between a pair (n1,n2) of overlapping reads
    """
    hang1 = n2["begin"] - n1["begin"]
    overlap = zip(n1["sites"][hang1:], n2["sites"])
    match, mismatch = (0, 0)
    for (c1, c2) in overlap:
        if c1 in ["A", "C", "G", "T"] and c2 in ["A", "C", "G", "T"]:
            if c1 == c2:
                match += 1
            else:
                mismatch += 1
    return (match, mismatch)
