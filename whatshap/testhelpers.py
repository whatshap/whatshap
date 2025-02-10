"""
Utility functions only used by unit tests
"""

import math
import textwrap
from collections import defaultdict
from whatshap.core import PhredGenotypeLikelihoods, Read, ReadSet, Genotype


def likelihoods_equal(a: PhredGenotypeLikelihoods, b: PhredGenotypeLikelihoods):
    for gt in a.genotypes():
        if not math.isclose(a[gt], b[gt], abs_tol=1e-9):
            return False
    return True


def string_to_readset(s, w=None, sample_ids=None, source_id=0, scale_quality=None):
    s = textwrap.dedent(s).strip()
    if w is not None:
        w = textwrap.dedent(w).strip().split("\n")
    rs = ReadSet()
    for index, line in enumerate(s.split("\n")):
        if len(line) == 0:
            continue
        if sample_ids is None:
            read = Read(f"Read {index + 1}", 50, source_id)
        else:
            read = Read(f"Read {index + 1}", 50, source_id, sample_ids[index])
        for pos, c in enumerate(line):
            if c == " ":
                continue
            q = 1
            if w is not None:
                q = int(w[index][pos])
            if scale_quality is not None:
                read.add_variant(position=(pos + 1) * 10, allele=int(c), quality=q * scale_quality)
            else:
                read.add_variant(position=(pos + 1) * 10, allele=int(c), quality=q)
        assert len(read) > 1, "Reads covering less than two variants are not allowed"
        rs.add(read)
    print(rs)
    return rs


def string_to_readset_pedigree(s, w=None, scaling_quality=None):
    s = textwrap.dedent(s).strip()
    read_sources = []
    s2 = ""
    for line in s.split("\n"):
        if len(line) == 0:
            continue
        individual = ord(line[0]) - ord("A")
        assert 0 <= individual < 26
        read_sources.append(individual)
        s2 += line[1:] + "\n"
    rs = string_to_readset(s=s2, w=w, sample_ids=read_sources, scale_quality=scaling_quality)
    print("read_sources:", read_sources)
    return rs


def matrix_to_readset(lines):
    rs = ReadSet()
    index_tracker = 0
    for line in lines:
        s = line.split()
        assert len(s) % 2 == 1, "Not in matrix format."

        index = int(s[0])
        index_tracker += 1
        assert index == index_tracker, "Not in matrix format."

        read = Read(f"Read {index}", 50)
        for i in range(int(len(s) / 2)):
            offset = int(s[2 * i + 1])
            for pos, c in enumerate(s[2 * i + 2]):
                read.add_variant(position=(offset + pos) * 10, allele=int(c), quality=1)

        rs.add(read)

    print(rs)
    return rs


def flip_cost(variant, target_value):
    """Returns cost of flipping the given read variant to target_value."""
    if variant.allele == target_value:
        return 0
    else:
        return variant.quality


def is_ambiguous(assignments):
    sets = [set(), set()]
    for assignment in assignments:
        for s, allele in zip(sets, assignment):
            s.add(allele)
    return [len(s) > 1 for s in sets]


def column_cost(variants, possible_assignments):
    """Compute cost for one position and return the minimum cost assignment.
    Returns ('X','X') if minimum is not unique (i.e. a "tie")."""
    costs = []
    for allele1, allele2 in possible_assignments:
        cost1 = sum(flip_cost(v, allele1) for v in variants[0])
        cost2 = sum(flip_cost(v, allele2) for v in variants[1])
        costs.append(cost1 + cost2)
    l = [(cost, i) for i, cost in enumerate(costs)]
    l.sort()
    min_cost = l[0][0]
    best_assignment = list(possible_assignments[l[0][1]])
    # check for ties
    counts = defaultdict(int)
    for cost, index in l:
        counts[cost] += 1
    ties = counts[min_cost]
    ambiguous = is_ambiguous([possible_assignments[i] for cost, i in l[:ties]])
    for i in range(2):
        if ambiguous[i]:
            best_assignment[i] = 3
    return min_cost, best_assignment


def brute_force_phase(read_set, all_heterozygous):
    """Solves MEC by enumerating all possible bipartitions."""

    def print(*args):
        pass

    assert len(read_set) < 10, "Too many reads for brute force"
    positions = read_set.get_positions()
    if all_heterozygous:
        possible_assignments = [(0, 1), (1, 0)]
    else:
        possible_assignments = [(0, 0), (0, 1), (1, 0), (1, 1)]
    # bit i in "partition" encodes to which set read i belongs
    best_partition = None
    best_cost = None
    best_haplotypes = None
    solution_count = 0
    for partition in range(2 ** len(read_set)):
        print(f"Looking at partition {partition:0>{len(read_set)}b}")
        # compute cost induced by that partition
        cost = 0
        haplotypes = []
        for p in positions:
            # find variants covering this position
            variants = [[], []]
            for n, read in enumerate(read_set):
                i = (partition >> n) & 1
                for variant in read:
                    if variant.position == p:
                        variants[i].append(variant)
            c, assignment = column_cost(variants, possible_assignments)
            print(f"    position: {p}, variants: {str(variants)} --> cost = {c}")
            cost += c
            haplotypes.append(assignment)
        print("  --> cost for this partitioning:", cost)
        if (best_cost is None) or (cost < best_cost):
            best_partition = partition
            best_cost = cost
            best_haplotypes = haplotypes
            solution_count = 1
        elif cost == best_cost:
            solution_count += 1
    # Each partition has its inverse with the same cost
    assert solution_count % 2 == 0
    haplotype1 = "".join([str(allele1) for allele1, allele2 in best_haplotypes])
    haplotype2 = "".join([str(allele2) for allele1, allele2 in best_haplotypes])
    return (
        best_cost,
        [(best_partition >> x) & 1 for x in range(len(read_set))],
        solution_count // 2,
        haplotype1,
        haplotype2,
    )


def canonic_index_to_biallelic_gt(num_alt, ploidy=2):
    """Takes the numeric VCF representation of a biallelic genotpyte and given ploidy
    Diploid:
    0 -> 0/0
    1 -> 0/1
    2 -> 1/1
    Trioploid:
    0 -> 0/0/0
    1 -> 0/0/1
    2 -> 0/1/1
    3 -> 1/1/1
    ...
    and converts it into a Genotype object

    See this link for further explanation:
    https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
    """
    if 0 <= num_alt <= ploidy:
        return Genotype([0] * (ploidy - num_alt) + [1] * (num_alt))
    else:
        return Genotype([])


def canonic_index_list_to_biallelic_gt_list(list_int, ploidy=2):
    """Returns a list of diploid, biallelic genotype objects
    according to the provided integer representation

    See this link for further explanation:
    https://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
    """
    return [canonic_index_to_biallelic_gt(i, ploidy) for i in list_int]
