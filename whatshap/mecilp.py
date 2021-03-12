from pulp import LpProblem, LpMinimize, LpVariable, LpContinuous, LpInteger, value


from whatshap.core import (
    Read,
    ReadSet,
)


def solve_mec(all_reads, accessible_positions, algorithm):
    index = {}
    rev_index = []
    num_vars = 0
    for position in all_reads.get_positions():
        index[position] = num_vars
        rev_index.append(position)
        num_vars += 1

    readlist = [dict() for _ in range(len(all_reads))]
    for i, read in enumerate(all_reads):
        for variant in read:
            readlist[i][index[variant.position]] = variant.allele

    reduce = True
    if algorithm == "ilpcompact":
        haplotypes, optimal_cost = mec_ilp_allhet_compact(readlist, reduce)
    elif algorithm == "ilpfast":
        haplotypes, optimal_cost = mec_ilp_allhet_fast(readlist, reduce)
    elif algorithm == "ilpfast-general":
        haplotypes, optimal_cost = mec_ilp_general_fast(readlist, reduce)
    superread = ReadSet()
    for i in range(2):
        read = Read("superread {}".format(i + 1), 0, 0)
        # insert alleles
        for j, allele in enumerate(haplotypes[i]):
            if allele is None:
                continue
            read.add_variant(accessible_positions[j], allele, 0)
        superread.add(read)
    superreads_list = [superread]
    transmission_vector = [0 for _ in range(num_vars)]

    return superreads_list, transmission_vector


def mec_ilp_allhet_compact(M):
    # M is a list of dictionaries, M[i][j] is allele of read i at position j

    # determine bounds
    n = len(M)
    m = 0
    for read in M:
        for j in read:
            m = max(m, j + 1)

    # compute zero-set and one-set of m for every column j
    zero_set = [[] for _ in range(m)]
    one_set = [[] for _ in range(m)]
    for i, read in enumerate(M):
        for pos in read:
            if read[pos] == 0:
                zero_set[pos].append(i)
            else:
                one_set[pos].append(i)

    # setup model
    model = LpProblem("MEC_allhet", LpMinimize)

    # constants
    c = [len(zero_set[j]) + len(one_set[j]) for j in range(m)]  # position-wise coverage

    # create variables #LpContinuous #LpInteger
    x = [None for _ in range(n)]
    for i in range(n):
        x[i] = LpVariable("x_{}".format(i), 0, 1, LpContinuous)  # haplotype to which read i goes
    h = [None for _ in range(m)]
    for j in range(m):
        h[j] = LpVariable("h_{}".format(j), 0, 1, LpInteger)  # allele of haplotype 0
    e = [None for _ in range(m)]
    for j in range(m):
        e[j] = LpVariable("e_{}".format(j), 0, None, LpContinuous)  # number of errors in column j

    # objective
    model += sum([e[j] for j in range(m)])

    # constraints
    for j in range(m):
        model += (
            e[j]
            >= sum([1 - x[i] for i in one_set[j]]) + sum([x[i] for i in zero_set[j]]) - c[j] * h[j]
        )
        model += e[j] >= sum([1 - x[i] for i in zero_set[j]]) + sum([x[i] for i in one_set[j]]) - c[
            j
        ] * (1 - h[j])
        # model += c[j]*h[j] <= c[j] - sum([1-x[i] for i in zero_set[j]]) - sum([x[i] for i in one_set[j]]) + sum([1-x[i] for i in one_set[j]]) + sum([x[i] for i in zero_set[j]])
        # model += c[j]*h[j] >= sum([1-x[i] for i in one_set[j]]) - sum([x[i] for i in zero_set[j]]) - sum([1-x[i] for i in zero_set[j]]) - sum([x[i] for i in one_set[j]])

    # reductions
    """
    if reduce:
        attached = 0
        for i in range(n):
            for j in range(i + 1, n):
                if covered[i][-1] <= covered[j][0]:
                    continue
                equal = 0
                diff = 0
                union_pos = set(covered[i]).union(set(covered[j]))
                inter_pos = set(covered[i]).intersection(set(covered[j]))
                excl = len(union_pos) - len(inter_pos)
                for p in inter_pos:
                    if M[i][p] == M[j][p]:
                        equal += 1
                    else:
                        diff += 1
                if equal > excl + diff:
                    model += x[i] == x[j]
                    attached += 1
        print("{} variable pairs attached.".format(attached))
    """

    model.solve()
    # model.solve(solvers.CPLEX_CMD(path="/opt/ibm/ILOG/CPLEX_Studio129/cplex/bin/x86-64_linux/cplex",msg=0))
    for xi in x:
        if xi.varValue is None or 0.001 < xi.varValue < 0.999:
            print("Warning: {} had non-integer assignment of {}.".format(xi, xi.varValue))
    for hj in h:
        if hj.varValue is None or 0.001 < hj.varValue < 0.999:
            print("Warning: {} had non-integer assignment of {}.".format(hj, hj.varValue))

    return variables_to_values(h), value(model.objective)


def mec_ilp_allhet_fast(M, reduce=False):
    # M is a list of dictionaries, M[i][j] is allele of read i at position j

    # determine bounds
    n = len(M)
    m = 0
    for read in M:
        for j in read:
            m = max(m, j + 1)

    # compute zero-set and one-set of m for every column j
    zero_set = [[] for _ in range(m)]
    one_set = [[] for _ in range(m)]
    covered = [[] for _ in range(n)]
    for i, read in enumerate(M):
        for pos in read:
            covered[i].append(pos)
            if read[pos] == 0:
                zero_set[pos].append(i)
            else:
                one_set[pos].append(i)

    # setup model
    model = LpProblem("MEC_allhet", LpMinimize)

    # create variables #LpContinuous #LpInteger
    x = [None for _ in range(n)]
    for i in range(n):
        x[i] = LpVariable("x_{}".format(i), 0, 1, LpContinuous)  # 1 iff read i goes to haplotype 0
    h = [None for _ in range(m)]
    for j in range(m):
        h[j] = LpVariable("h_{}".format(j), 0, 1, LpInteger)  # allele of haplotype 0
    f = [{} for i in range(n)]
    for i in range(n):
        for j in covered[i]:
            f[i][j] = LpVariable(
                "f_{}_{}".format(i, j), 0, None, LpContinuous
            )  # 1 iff read i goes to hap 0 but is different from it

    # objective
    model += sum(
        [sum([(1 - h[j] - x[i] + 2 * f[i][j]) for i in zero_set[j]]) for j in range(m)]
    ) + sum([sum([(h[j] - x[i] + 2 * f[i][j]) for i in one_set[j]]) for j in range(m)])

    # constraints
    for j in range(m):
        for i in zero_set[j]:
            model += h[j] + x[i] - 1 <= f[i][j]
        for i in one_set[j]:
            model += x[i] - h[j] <= f[i][j]

    # reductions
    if reduce:
        attached = 0
        for i in range(n):
            for j in range(i + 1, n):
                if covered[i][-1] <= covered[j][0]:
                    continue
                equal = 0
                diff = 0
                union_pos = set(covered[i]).union(set(covered[j]))
                inter_pos = set(covered[i]).intersection(set(covered[j]))
                excl = len(union_pos) - len(inter_pos)
                for p in inter_pos:
                    if M[i][p] == M[j][p]:
                        equal += 1
                    else:
                        diff += 1
                if equal > excl + diff:
                    model += x[i] == x[j]
                    attached += 1
        print("{} variable pairs attached.".format(attached))

    model.solve()
    for xi in x:
        if 0.001 < xi.varValue < 0.999:
            print("Warning: {} had non-integer assignment of {}.".format(xi, xi.varValue))
    for hj in h:
        if hj.varValue is None or 0.001 < hj.varValue < 0.999:
            print("Warning: {} had non-integer assignment of {}.".format(hj, hj.varValue))
    for i in range(n):
        for j in covered[i]:
            if f[i][j].varValue is None or 0.001 < f[i][j].varValue < 0.999:
                print(
                    "Warning: {} had non-integer assignment of {}.".format(
                        f[i][j], f[i][j].varValue
                    )
                )

    return variables_to_values(h), value(model.objective)


def mec_ilp_general_fast(M, reduce=False):
    # M is a list of dictionaries, M[i][j] is allele of read i at position j

    # determine bounds
    n = len(M)
    m = 0
    for read in M:
        for j in read:
            m = max(m, j + 1)

    # compute zero-set and one-set of m for every column j
    zero_set = [[] for _ in range(m)]
    one_set = [[] for _ in range(m)]
    covered = [[] for _ in range(n)]
    for i, read in enumerate(M):
        for pos in read:
            covered[i].append(pos)
            if read[pos] == 0:
                zero_set[pos].append(i)
            else:
                one_set[pos].append(i)

    # setup model
    model = LpProblem("MEC_allhet", LpMinimize)

    # create variables #LpContinuous #LpInteger
    x = [None for _ in range(n)]
    for i in range(n):
        x[i] = LpVariable("x_{}".format(i), 0, 1, LpContinuous)  # 1 iff read i goes to haplotype 0
    h0 = [None for _ in range(m)]
    h1 = [None for _ in range(m)]
    for j in range(m):
        h0[j] = LpVariable("h0_{}".format(j), 0, 1, LpInteger)  # allele of haplotype 0
        h1[j] = LpVariable("h1_{}".format(j), 0, 1, LpInteger)  # allele of haplotype 0
    f = [{} for i in range(n)]
    for i in range(n):
        for j in covered[i]:
            f[i][j] = LpVariable(
                "f_{}_{}".format(i, j), 0, None, LpContinuous
            )  # 1 iff read i goes to hap 0 but is different from it

    # objective
    model += sum([sum([f[i][j] for j in covered[i]]) for i in range(n)])
    # model += sum([sum([f[i][j] for i in zero_set[j]]) for j in range(m)]) + sum([sum([f[i][j] for i in one_set[j]]) for j in range(m)])

    # constraints
    for j in range(m):
        for i in zero_set[j]:
            model += h0[j] + x[i] - 1 <= f[i][j]
            model += h1[j] - x[i] <= f[i][j]
        for i in one_set[j]:
            model += x[i] - h0[j] <= f[i][j]
            model += 1 - x[i] - h1[j] <= f[i][j]

    # reductions
    if reduce:
        attached = 0
        for i in range(n):
            for j in range(i + 1, n):
                if covered[i][-1] <= covered[j][0]:
                    continue
                equal = 0
                diff = 0
                union_pos = set(covered[i]).union(set(covered[j]))
                inter_pos = set(covered[i]).intersection(set(covered[j]))
                excl = len(union_pos) - len(inter_pos)
                for p in inter_pos:
                    if M[i][p] == M[j][p]:
                        equal += 1
                    else:
                        diff += 1
                if equal > excl + diff:
                    model += x[i] == x[j]
                    attached += 1
        print("{} variable pairs attached.".format(attached))

    model.solve()
    for xi in x:
        if 0.001 < xi.varValue < 0.999:
            print("Warning: {} had non-integer assignment of {}.".format(xi, xi.varValue))
    for hj in h0:
        if hj.varValue is None or 0.001 < hj.varValue < 0.999:
            print("Warning: {} had non-integer assignment of {}.".format(hj, hj.varValue))
    for hj in h1:
        if hj.varValue is None or 0.001 < hj.varValue < 0.999:
            print("Warning: {} had non-integer assignment of {}.".format(hj, hj.varValue))
    for i in range(n):
        for j in covered[i]:
            if f[i][j].varValue is None or 0.001 < f[i][j].varValue < 0.999:
                print(
                    "Warning: {} had non-integer assignment of {}.".format(
                        f[i][j], f[i][j].varValue
                    )
                )

    haps = variables_to_values_general(h0, h1)
    for j in range(m):
        if haps[0][j] == haps[1][j]:
            print("Found homozygous site: {}".format(j))
            print(zero_set[j])
            print(one_set[j])

    return haps, value(model.objective)


def variables_to_values(hap_var_list):
    m = len(hap_var_list)
    hap0 = []
    hap1 = []
    for j in range(m):
        if hap_var_list[j].varValue is None:
            hap0.append(None)
            hap1.append(None)
        else:
            if hap_var_list[j].varValue > 0.999:
                hap0.append(1)
                hap1.append(0)
            else:
                hap0.append(0)
                hap1.append(1)

    return [hap0, hap1]


def variables_to_values_general(h0, h1):
    m = len(h0)
    hap0 = []
    hap1 = []
    for j in range(m):
        if h0[j].varValue is None:
            hap0.append(None)
        else:
            if h0[j].varValue > 0.999:
                hap0.append(1)
            else:
                hap0.append(0)
        if h1[j].varValue is None:
            hap1.append(None)
        else:
            if h1[j].varValue > 0.999:
                hap1.append(1)
            else:
                hap1.append(0)

    return [hap0, hap1]
