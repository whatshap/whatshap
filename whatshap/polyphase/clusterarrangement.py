"""
For a clustering of marker alleles in the context of polyploid genetic phasing, this module
computes an assignment of these alleles to a set of haplotypes in order to obtain a (partial)
phasing of the variants to which these alleles belong.
"""

import logging

from pulp import LpProblem, LpVariable, LpMaximize, LpInteger, value
from whatshap.polyphase import get_ilp_solver

logger = logging.getLogger(__name__)


def arrange_clusters(clustering, padding, ploidy):
    # filter out singleton clusters
    filtered_to_real = []
    fclustering = []
    for i, clust in enumerate(clustering):
        if len(clust) > 1:
            filtered_to_real.append(i)
            fclustering.append(clust)

    # determine start and end of clusters (including padding) and worth (=number of covered variants)
    c_start = []
    c_end = []
    c_worth = []
    for clust in fclustering:
        variants = [v for v in clust]
        c_worth.append(len(variants))
        c_start.append(max(0, min(variants) - padding))
        c_end.append(max(variants) + padding)

    n = max(c_end)
    c = len(fclustering)

    # setup model
    model = LpProblem("Cluster_Arrangement_c{}_n{}_p{}".format(c, n, ploidy), LpMaximize)

    # x[i][j] = 1 if cluster j is put on haplotype i, else 0
    x = [
        [LpVariable("x_{}_{}".format(i, j), 0, 1, LpInteger) for j in range(c)]
        for i in range(ploidy)
    ]

    # maximize worth of picked clusters (= maximize number of explained variants)
    model += sum([c_worth[j] * x[i][j] for j in range(c) for i in range(ploidy)])

    # cluster can only be put on at most one haplotype
    for j in range(c):
        model += sum([x[i][j] for i in range(ploidy)]) <= 1

    # for each position: at most one overlapping cluster can be on same haplotype
    old_covered = []
    for pos in range(n):
        covered = sorted([i for i in range(c) if c_start[i] <= pos <= c_end[i]])
        if covered != old_covered:
            for i in range(ploidy):
                model += sum([x[i][j] for j in covered]) <= 1
            old_covered = covered

    # solve model
    solver = get_ilp_solver()
    model.solve(solver)

    selected = []

    objVal = value(model.objective)
    logger.info(
        "Arranged %i variants out of a total of %i",
        int(objVal),
        sum([len(clust) for clust in clustering]),
    )

    for i in range(ploidy):
        selected.append([filtered_to_real[j] for j in range(c) if x[i][j].varValue > 0.999])
        logger.info("   h%i: %s", i, selected[-1])

    return selected
