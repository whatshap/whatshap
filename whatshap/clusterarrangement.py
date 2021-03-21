import logging

from pulp import LpProblem, LpVariable, LpMaximize, LpContinuous, LpInteger, value, COIN_CMD

#from whatshap.vcf import VariantTable

logger = logging.getLogger(__name__)


def arrange_clusters(clustering, node_to_variant, padding, ploidy):

    # determine start and end of clusters (including padding) and worth (=number of covered variants)
    c_start = []
    c_end = []
    c_worth = []
    for clust in clustering:
        #variants = [node_to_variant[v] for v in clust]
        variants = [v for v in clust]
        c_worth.append(len(variants))
        c_start.append(max(0, min(variants) - padding))
        c_end.append(max(variants) + padding)
        
    n = max(c_end)
    c = len(clustering)
        
    # setup model
    model = LpProblem("Cluster_Arrangement_c{}_n{}_p{}".format(c, n, ploidy), LpMaximize)
    
    # x[i][j] = 1 if cluster j is put on haplotype i, else 0
    x = [[LpVariable("x_{}_{}".format(i, j), 0, 1, LpInteger) for j in range(c)] for i in range(ploidy)]
    
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
    solver = COIN_CMD(msg=0)
    model.solve()
    
    selected = []
    
    objVal = value(model.objective)
    print("Arranged {} variants out of a total of {}".format(objVal, sum([len(clust) for clust in clustering])))
    
    for i in range(ploidy):
        selected.append([j for j in range(c) if x[i][j].varValue > 0.999])
        print("h{}: {}".format(i, selected[-1]))
    
    return selected
