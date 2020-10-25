import logging

from math import log
from typing import List
from scipy.stats import binom

from whatshap.core import (
    Genotype,
    TriangleSparseMatrix,
)
from whatshap.polyphaseplots import create_histogram
from whatshap.vcf import VariantTable

logger = logging.getLogger(__name__)


def get_variant_scoring(variant_table: VariantTable, parent: str, co_parent: str, offspring: List[str]):
    scoring = TriangleSparseMatrix()
    
    num_vars = len(variant_table)
    max_dist = 80
    
    ref, alt, alt_count, alt_count_co = classify_variants(variant_table, parent, co_parent)
    
    num_nodes = 0
    node_to_variant = dict()
    node_positions = []
    simplex_nulliplex_nodes = []
    allowed_pairs = [(1,0)]
    for i in range(len(variant_table)):
        if alt[i] is None:
            continue
        if alt_count[i] == 1 and alt_count_co[i] == 0:
            simplex_nulliplex_nodes.append(num_nodes)
        if (alt_count[i], alt_count_co[i]) not in allowed_pairs:
            continue
        for j in range(alt_count[i]):
            node_to_variant[num_nodes] = i
            node_positions.append(i)
            num_nodes += 1
    
    logger.info("Number of nodes to cluster: %d", num_nodes)
    logger.info("Number of simplex-nulliplex variants: %d", len(simplex_nulliplex_nodes))
    
    frac_0_0 = []
    
    for i in range(num_nodes):
        ni = node_to_variant[i]
        # skip if position i is not simplex-nulliplex
        if alt_count[ni] != 1 or alt_count_co[ni] != 0:
            continue
        
        print("scoring node {}: ref={}, alt={}, count = {} / {}".format(i,ref[ni], alt[ni], alt_count[ni], alt_count_co[ni]))
        # iterate over next max_dist relevant positions
        for j in range(i+1, min(i+max_dist+1, num_nodes)):
            nj = node_to_variant[j]
            if node_to_variant[i] == node_to_variant[j]:
                score = -float("inf")
                print("    partner {} on same variant with score".format(j, score))
            else:
                off_gts = [(variant_table.genotypes_of(off)[ni], variant_table.genotypes_of(off)[nj]) for off in offspring]
                off_gts = []
                for off in offspring:
                    gt = variant_table.genotypes_of(off)
                    if gt[ni].is_none() or gt[nj].is_none():
                        continue
                    else:
                        off_gts.append((gt[ni], gt[nj]))
                
                score = 0.0
                
                if len(off_gts) > 0:
                    if alt_count[nj] == 1 and alt_count_co[nj] == 0:
                        score = score_simplex_nulliplex_tetra(off_gts, ref_1=ref[ni], alt_1=alt[ni],  ref_2=ref[nj], alt_2=alt[nj])
                        
                        target_gts = (Genotype([ref[ni]]*4), Genotype([ref[nj]]*4))
                        count_0_0 = count_target_gt(off_gts, target_gts)
                        frac_0_0.append(count_0_0/len(off_gts))
                        
                    elif alt_count[nj] == 2 and alt_count_co[nj] == 0:
                        score = score_duplex_nulliplex_tetra(off_gts, ref_1=ref[ni], alt_1=alt[ni],  ref_2=ref[nj], alt_2=alt[nj])
                    elif alt_count[nj] == 3 and alt_count_co[nj] == 0:
                        score = score_triplex_nulliplex_tetra(off_gts, ref_1=ref[ni], alt_1=alt[ni],  ref_2=ref[nj], alt_2=alt[nj])
                    elif alt_count[nj] == 1 and alt_count_co[nj] == 1:
                        score = score_simplex_simplex_tetra(off_gts, ref_1=ref[ni], alt_1=alt[ni],  ref_2=ref[nj], alt_2=alt[nj])
                '''
                if score > 0.0:
                    score = 1.0
                elif score < 0.0:
                    score = -1.0
                '''
                print("    partner {}: ref={}, alt={}, count = {} / {} -> {} (using {} children)".format(j, ref[nj], alt[nj], alt_count[nj], alt_count_co[nj], score, len(off_gts)))
            
            scoring.set(i, j, score)
    
    create_histogram("fractions.pdf", frac_0_0, [], 50, [0.0, 1.0], "fraction", "Fraction of 0/0/0/0-0/0/0/0 in Simplex-Nulliplex", name1='0-0', name2='n/a')
            
    return scoring, node_to_variant


def classify_variants(variant_table: VariantTable, parent: str, co_parent: str):
    
    ref, alt, alt_count, alt_count_co = dict(), dict(), dict(), dict()
    gts1 = variant_table.genotypes_of(parent)
    gts2 = variant_table.genotypes_of(co_parent)
    
    for i in range(len(variant_table)):
        gt1 = gts1[i]
        gt2 = gts2[i]
        gt1v = gt1.as_vector()
        gt2v = gt2.as_vector()
        
        if gt1.is_none() or gt2.is_none():
            ref[i] = None
            alt[i] = None
            alt_count[i] = 0
            alt_count_co[i] = 0
            continue
        
        # check homozygosity
        if gt1.is_homozygous():
            ref[i] = gt1v[0]
            alt[i] = None
            alt_count[i] = 0
            alt_count_co[i] = 0
            continue

        # count different alleles
        alleles = set()
        for gt in [gt1v, gt2v]:
            for a in gt:
                alleles.add(a)

        alleles = sorted(list(alleles))

        if len(alleles) > 2:
            # genotypes are not bi-allelic
            ref[i] = None
            alt[i] = None
            alt_count[i] = 0
            alt_count_co[i] = 0
            continue

        assert len(alleles) == 2
        
        # determine majority allele and its multiplicity
        gt1v.sort()
        ref[i] = gt1v[int(len(gt1v)/2-1)]
        alt[i] = gt1v[0] if gt1v[0] != ref[i] else gt1v[-1]
        alt_count[i] = sum([1 if a == alt[i] else 0 for a in gt1v])
        alt_count_co[i] = sum([1 if a == alt[i] else 0 for a in gt2v])
        
    return ref, alt, alt_count, alt_count_co
    

def score_simplex_nulliplex_tetra(off_gts, ref_1=0, alt_1=1, ref_2=0, alt_2=1):
    target_gts = (Genotype([ref_1]*4), Genotype([ref_2]*4))
    count_0_0 = count_target_gt(off_gts, target_gts)
    
    # if both variants share 1-allele: support = 1/2
    pval_e = binom.pmf(count_0_0, len(off_gts), 1/2)
    
    # if both variants don't share 1-allele: support = 1/6
    pval_d = binom.pmf(count_0_0, len(off_gts), 1/6)
    
    print("    eval: {} out of {}: log({}/{}) = {}".format(count_0_0, len(off_gts), pval_e, pval_d, log(pval_e/pval_d)))
    
    return log(pval_e / pval_d)


def score_duplex_nulliplex_tetra(off_gts, ref_1=0, alt_1=1, ref_2=0, alt_2=1):
    target_gts = (Genotype([ref_1]*4), Genotype([ref_2]*4))
    count_0_0 = count_target_gt(off_gts, target_gts)
    
    # if both variants share 1-allele: support = 1/6
    pval_e = binom.pmf(count_0_0, len(off_gts), 1/6)
    
    # if both variants don't share 1-allele: support = 0
    pval_d = binom.pmf(count_0_0, len(off_gts), 1/24)
    
    return log(pval_e / pval_d)


def score_triplex_nulliplex_tetra(off_gts, ref_1=0, alt_1=1, ref_2=0, alt_2=1):
    target_gts = (Genotype([alt_1]*4), Genotype([alt_2]*4))
    count_0_0 = count_target_gt(off_gts, target_gts)
    
    # if both variants share 1-allele: support = 1/6
    pval_e = binom.pmf(count_0_0, len(off_gts), 1/6)
    
    # if both variants don't share 1-allele: support = 1/2
    pval_d = binom.pmf(count_0_0, len(off_gts), 1/2)
    
    return log(pval_e / pval_d)


def score_simplex_simplex_tetra(off_gts, ref_1=0, alt_1=1, ref_2=0, alt_2=1):
    target_gts = (Genotype([ref_1]*4), Genotype([ref_2]*4))
    count_0_0 = count_target_gt(off_gts, target_gts)
    
    # if both variants share 1-allele: support = 1/4
    pval_e = binom.pmf(count_0_0, len(off_gts), 1/4)
    
    # if both variants don't share 1-allele: support = 1/12
    pval_d = binom.pmf(count_0_0, len(off_gts), 1/12)
    
    return log(pval_e / pval_d)


def count_target_gt(off_gts, target_gts):
    count = 0
    for gt1, gt2 in off_gts:
        if gt1 == target_gts[0] and gt2 == target_gts[1]:
            count += 1
    return count