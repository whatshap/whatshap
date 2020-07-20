"""
Test Threading
"""

from collections import defaultdict
from whatshap.threading import (
    get_position_map,
    get_coverage,
    get_coverage_absolute,
    get_cluster_start_end_positions,
    compute_cut_positions,
    improve_path_on_multiswitches,
    compute_threading_path,
    get_pos_to_clusters_map,
    get_local_cluster_consensus,
)
from whatshap.core import Read, ReadSet


def create_testinstance1():
    var_pos = [
        24,
        56,
        89,
        113,
        162,
        166,
        187,
        205,
        211,
        248,
        273,
        299,
        307,
        324,
        351,
        370,
        378,
        400,
        441,
        455,
        478,
        492,
    ]

    readset = ReadSet()

    matrix = [
        "0011000",
        "11010100",
        " 101011010",
        " 0001011000",
        "  11001001",
        "  0010100000",
        "   100010001",
        "       0100000101",
        "    101110001",
        "        0001110011",
        "        1010001010",
        "     011100011",
        "         0010100111",
        "          1010101011",
        "          0101001110",
        "              01000001",
        "              01010001",
        "                101100",
        "                111010",
    ]

    for i in range(len(matrix)):
        read = Read(name="read" + str(i), mapq=15)
        for j in range(len(matrix[i])):
            if matrix[i][j] != " ":
                read.add_variant(var_pos[j], int(matrix[i][j]), 0)
        readset.add(read)

    clustering = [
        [0, 4, 6],
        [1, 2],
        [7, 10, 13],
        [9, 12, 14],
        [3, 5, 8, 11],
        [15, 16],
        [17],
        [18],
    ]
    genotypes = [
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 1, 1: 2},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 3, 1: 0},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 1, 1: 2},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 1, 1: 2},
        {0: 2, 1: 1},
        {0: 1, 1: 2},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
        {0: 2, 1: 1},
    ]

    return readset, var_pos, clustering, genotypes


def test_auxiliary_datastructures():
    # test postion map
    readset, var_pos, _, _ = create_testinstance1()
    index, rev_index = get_position_map(readset)
    for i in range(len(var_pos)):
        assert index[var_pos[i]] == i
    assert rev_index == var_pos

    # test relative coverage
    clustering = [
        [0, 4, 6],
        [1, 2],
        [7, 10, 13],
        [9, 12, 14],
        [3, 5, 8, 11],
        [15, 16],
        [17],
        [18],
    ]
    cov = get_coverage(readset, clustering, index)
    assert cov[0] == {0: 0.5, 1: 0.5}
    assert cov[1] == {0: 0.25, 1: 0.5, 4: 0.25}
    assert cov[2] == {0: 1 / 3, 1: 1 / 3, 4: 1 / 3}
    assert cov[3] == {0: 3 / 7, 1: 2 / 7, 4: 2 / 7}
    assert cov[4] == {0: 3 / 8, 1: 2 / 8, 4: 3 / 8}
    assert cov[5] == {0: 3 / 9, 1: 2 / 9, 4: 4 / 9}
    assert cov[6] == {0: 3 / 9, 1: 2 / 9, 4: 4 / 9}
    assert cov[7] == {0: 2 / 9, 1: 2 / 9, 2: 1 / 9, 4: 4 / 9}
    assert cov[8] == {0: 2 / 10, 1: 1 / 10, 2: 2 / 10, 3: 1 / 10, 4: 4 / 10}
    assert cov[9] == {0: 2 / 11, 1: 1 / 11, 2: 2 / 11, 3: 2 / 11, 4: 4 / 11}
    assert cov[10] == {0: 1 / 11, 2: 3 / 11, 3: 3 / 11, 4: 4 / 11}
    assert cov[11] == {0: 1 / 10, 2: 3 / 10, 3: 3 / 10, 4: 3 / 10}
    assert cov[12] == {2: 3 / 8, 3: 3 / 8, 4: 2 / 8}
    assert cov[13] == {2: 3 / 7, 3: 3 / 7, 4: 1 / 7}
    assert cov[14] == {2: 3 / 8, 3: 3 / 8, 5: 2 / 8}
    assert cov[15] == {2: 3 / 8, 3: 3 / 8, 5: 2 / 8}
    assert cov[16] == {2: 3 / 10, 3: 3 / 10, 5: 2 / 10, 6: 1 / 10, 7: 1 / 10}
    assert cov[17] == {2: 2 / 9, 3: 3 / 9, 5: 2 / 9, 6: 1 / 9, 7: 1 / 9}
    assert cov[18] == {2: 1 / 7, 3: 2 / 7, 5: 2 / 7, 6: 1 / 7, 7: 1 / 7}
    assert cov[19] == {2: 1 / 6, 3: 1 / 6, 5: 2 / 6, 6: 1 / 6, 7: 1 / 6}
    assert cov[20] == {5: 2 / 4, 6: 1 / 4, 7: 1 / 4}
    assert cov[21] == {5: 2 / 4, 6: 1 / 4, 7: 1 / 4}

    # test absolute coverage
    abs_cov = get_coverage_absolute(readset, clustering, index)
    assert abs_cov[0] == {0: 1, 1: 1}
    assert abs_cov[1] == {0: 1, 1: 2, 4: 1}
    assert abs_cov[2] == {0: 2, 1: 2, 4: 2}
    assert abs_cov[3] == {0: 3, 1: 2, 4: 2}
    assert abs_cov[4] == {0: 3, 1: 2, 4: 3}
    assert abs_cov[5] == {0: 3, 1: 2, 4: 4}
    assert abs_cov[6] == {0: 3, 1: 2, 4: 4}
    assert abs_cov[7] == {0: 2, 1: 2, 2: 1, 4: 4}
    assert abs_cov[8] == {0: 2, 1: 1, 2: 2, 3: 1, 4: 4}
    assert abs_cov[9] == {0: 2, 1: 1, 2: 2, 3: 2, 4: 4}
    assert abs_cov[10] == {0: 1, 2: 3, 3: 3, 4: 4}
    assert abs_cov[11] == {0: 1, 2: 3, 3: 3, 4: 3}
    assert abs_cov[12] == {2: 3, 3: 3, 4: 2}
    assert abs_cov[13] == {2: 3, 3: 3, 4: 1}
    assert abs_cov[14] == {2: 3, 3: 3, 5: 2}
    assert abs_cov[15] == {2: 3, 3: 3, 5: 2}
    assert abs_cov[16] == {2: 3, 3: 3, 5: 2, 6: 1, 7: 1}
    assert abs_cov[17] == {2: 2, 3: 3, 5: 2, 6: 1, 7: 1}
    assert abs_cov[18] == {2: 1, 3: 2, 5: 2, 6: 1, 7: 1}
    assert abs_cov[19] == {2: 1, 3: 1, 5: 2, 6: 1, 7: 1}
    assert abs_cov[20] == {5: 2, 6: 1, 7: 1}
    assert abs_cov[21] == {5: 2, 6: 1, 7: 1}


def test_clusterbased_structures():
    readset, var_pos, clustering, _ = create_testinstance1()
    index, rev_index = get_position_map(readset)

    # clustering bounds
    cluster_start_ends = get_cluster_start_end_positions(readset, clustering, index)
    assert cluster_start_ends[0] == (0, 11)
    assert cluster_start_ends[1] == (0, 9)
    assert cluster_start_ends[2] == (7, 19)
    assert cluster_start_ends[3] == (8, 19)
    assert cluster_start_ends[4] == (1, 13)
    assert cluster_start_ends[5] == (14, 21)
    assert cluster_start_ends[6] == (16, 21)
    assert cluster_start_ends[7] == (16, 21)


def test_cut_positions():
    path = [
        [2, 3, 5, 1],
        [2, 3, 5, 1],
        [2, 3, 6, 1],
        [2, 3, 6, 1],
        [4, 8, 6, 1],
        [5, 8, 6, 1],
        [5, 8, 6, 1],
        [5, 8, 6, 6],
        [8, 8, 6, 6],
        [8, 8, 6, 7],
        [8, 8, 6, 7],
        [8, 8, 9, 10],
        [8, 11, 9, 10],
        [8, 11, 9, 10],
    ]

    cuts1 = compute_cut_positions(path, 1, 12)
    cuts2 = compute_cut_positions(path, 2, 12)
    cuts3 = compute_cut_positions(path, 3, 12)
    cuts4 = compute_cut_positions(path, 4, 12)
    cuts5 = compute_cut_positions(path, 5, 12)

    assert cuts1[0] == [0]
    assert cuts2[0] == [0]
    assert cuts3[0] == [0, 4, 11]
    assert cuts4[0] == [0, 4, 9, 11]
    assert cuts5[0] == [0, 2, 4, 5, 7, 8, 9, 11, 12]

    assert cuts1[1] == [[0], [0], [0], [0]]
    assert cuts2[1] == [[0], [0], [0], [0]]
    assert cuts3[1] == [[0, 4], [0, 4], [0, 11], [0, 11]]
    assert cuts4[1] == [[0, 4], [0, 4], [0, 9, 11], [0, 9, 11]]
    assert cuts5[1] == [[0, 4, 5, 8, 12], [0, 4, 12], [0, 2, 9, 11], [0, 7, 9, 11]]


def test_multiswitch_improvement():
    path = [
        [3, 1, 2, 4],
        [3, 1, 2, 4],
        [5, 1, 2, 4],
        [5, 1, 2, 4],
        [5, 7, 6, 4],
        [5, 7, 6, 4],
        [5, 7, 6, 7],
        [5, 7, 6, 4],
        [5, 7, 6, 4],
        [8, 9, 10, 4],
        [8, 9, 10, 4],
    ]
    cluster_sim = [defaultdict(float) for _ in range(len(path))]
    for i in range(len(path)):
        cluster_sim[i][(1, 7)] = 0.6
        cluster_sim[i][(1, 6)] = 0.7
        cluster_sim[i][(2, 7)] = 0.8
        cluster_sim[i][(2, 6)] = 0.65
        cluster_sim[i][(5, 8)] = 0.3
        cluster_sim[i][(5, 9)] = 0.5
        cluster_sim[i][(5, 10)] = 0.8
        cluster_sim[i][(7, 8)] = 0.5
        cluster_sim[i][(7, 9)] = 0.8
        cluster_sim[i][(7, 10)] = 0.85
        cluster_sim[i][(6, 8)] = 0.9
        cluster_sim[i][(6, 9)] = 0.9
        cluster_sim[i][(6, 10)] = 0.6

    corrected_path = improve_path_on_multiswitches(path, 11, cluster_sim)

    truth = [
        [3, 1, 2, 4],
        [3, 1, 2, 4],
        [5, 1, 2, 4],
        [5, 1, 2, 4],
        [5, 6, 7, 4],
        [5, 6, 7, 4],
        [5, 6, 7, 7],
        [5, 6, 7, 4],
        [5, 6, 7, 4],
        [10, 8, 9, 4],
        [10, 8, 9, 4],
    ]

    for i in range(len(truth)):
        assert corrected_path[i] == truth[i]


def test_path_no_affine():
    readset, var_pos, clustering, genotypes = create_testinstance1()
    ploidy = 3

    index, rev_index = get_position_map(readset)
    num_vars = len(rev_index)
    positions = get_cluster_start_end_positions(readset, clustering, index)
    coverage = get_coverage(readset, clustering, index)
    cov_map = get_pos_to_clusters_map(coverage, ploidy)
    consensus = get_local_cluster_consensus(readset, clustering, cov_map, positions)

    path = compute_threading_path(
        readset,
        clustering,
        num_vars,
        coverage,
        cov_map,
        consensus,
        ploidy,
        genotypes,
        affine_switch_cost=0.0,
    )
    cluster_paths = ["".join([str(path[i][j]) for i in range(len(path))]) for j in range(3)]

    first_block = set([cluster_paths[0][:20], cluster_paths[1][:20], cluster_paths[2][:20]])
    first_truth = set(["00000000003333333333", "11111111222222222222", "04444444444444555555"])
    second_block = set([cluster_paths[0][20:], cluster_paths[1][20:], cluster_paths[2][20:]])
    second_truth = set(["66", "77", "55"])

    print(cluster_paths)

    assert first_block == first_truth
    assert second_block == second_truth


def test_path_with_affine():
    readset, var_pos, clustering, genotypes = create_testinstance1()
    ploidy = 3

    index, rev_index = get_position_map(readset)
    num_vars = len(rev_index)
    positions = get_cluster_start_end_positions(readset, clustering, index)
    coverage = get_coverage(readset, clustering, index)
    cov_map = get_pos_to_clusters_map(coverage, ploidy)
    consensus = get_local_cluster_consensus(readset, clustering, cov_map, positions)

    path = compute_threading_path(
        readset, clustering, num_vars, coverage, cov_map, consensus, ploidy, genotypes
    )
    cluster_paths = ["".join([str(path[i][j]) for i in range(len(path))]) for j in range(3)]

    first_block = set([cluster_paths[0][:9], cluster_paths[1][:9], cluster_paths[2][:9]])
    first_truth = set(["000000000", "111111111", "044444444"])
    second_block = set([cluster_paths[0][9:20], cluster_paths[1][9:20], cluster_paths[2][9:20]])
    second_truth = set(["33333333333", "22222222222", "44444555555"])
    third_block = set([cluster_paths[0][20:], cluster_paths[1][20:], cluster_paths[2][20:]])
    third_truth = set(["66", "77", "55"])

    print(cluster_paths)

    assert first_block == first_truth
    assert second_block == second_truth
    assert third_block == third_truth
