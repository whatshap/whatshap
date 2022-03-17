# flake8: noqa
import itertools as it
from math import ceil, floor
from copy import deepcopy

import logging
from whatshap.core import Read, ReadSet
from whatshap.cli.compare import compute_switch_flips_poly_bt
from whatshap.polyphaseutil import get_coverage
from collections import defaultdict

"""
This class is exclusively used for debugging and development.
"""

logger = logging.getLogger(__name__)


def draw_plots(
    block_readsets,
    clustering,
    threading,
    haplotypes,
    cut_positions,
    phasable_variant_table,
    plot_clusters,
    plot_threading,
    output,
):
    # Plot options
    logger.info("Generating plots ...")
    combined_readset = ReadSet()
    for block_readset in block_readsets:
        for read in block_readset:
            combined_readset.add(read)
    if plot_clusters:
        draw_clustering(
            combined_readset,
            clustering,
            phasable_variant_table,
            output + ".clusters.pdf",
            genome_space=False,
        )
    if plot_threading:
        coverage = get_coverage(combined_readset, clustering)
        draw_threading(
            combined_readset,
            clustering,
            coverage,
            threading,
            cut_positions,
            haplotypes,
            phasable_variant_table,
            output + ".threading.pdf",
        )


"""
This method only works for a test dataset, for which the true haplotype of read was encoded
into its name. For any other read name, it just returns -1 for unknown haplotype
"""


def parse_haplotype(name):
    try:
        tokens = name.split("_")
        if tokens[-2] == "HG00514" and tokens[-1] == "HAP1":
            return 0
        elif tokens[-2] == "HG00514" and tokens[-1] == "HAP2":
            return 1
        elif tokens[-2] == "NA19240" and tokens[-1] == "HAP1":
            return 2
        elif tokens[-2] == "NA19240" and tokens[-1] == "HAP2":
            return 3
        elif tokens[-2] == "HG00733" and tokens[-1] == "HAP1":
            return 4
        elif tokens[-2] == "HG00733" and tokens[-1] == "HAP2":
            return 5
    except:
        pass
    return -1


def avg_readlength(readset):
    # Estiamtes the average read length in base pairs
    if len(readset) > 0:
        return sum([read[-1].position - read[0].position for read in readset]) / len(readset)
    else:
        return 0


def get_phase(readset, var_table):
    tmp_table = deepcopy(var_table)
    tmp_table.subset_rows_by_position(readset.get_positions())
    try:
        phase_rows = [variant.phase for variant in tmp_table.phases[0]]
    except AttributeError as e:
        return None
    return [[row[i] for row in phase_rows] for i in range(len(phase_rows[0]))]


def draw_clustering(readset, clustering, var_table, path, genome_space=False):
    try:
        import matplotlib

        matplotlib.use("agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from pylab import savefig

        # Sort a deep copy of clustering
        clusters = sorted(
            deepcopy(clustering),
            key=lambda x: min([readset[i][0].position for i in x]) if len(x) > 0 else 0,
        )

        # Construct real and predicted haplotype per read
        true_hap = [parse_haplotype(read.name) for read in readset]

        # Map variant positions to [0,l)
        index = {}
        rev_index = []
        num_vars = 0
        min_pos = float("inf")
        max_pos = 0
        for position in readset.get_positions():
            index[position] = num_vars
            rev_index.append(position)
            num_vars += 1
            min_pos = min(min_pos, position)
            max_pos = max(max_pos, position)

        min_pos = min(readset.get_positions()) if genome_space else 0
        max_pos = max(readset.get_positions()) if genome_space else num_vars

        # Plot heatmaps
        fig = plt.figure(figsize=(num_vars / 40, len(readset) / 40), dpi=100)
        legend_handles = {}
        y_offset = 0
        y_margin = 5

        # Plot haplotype dissimilarity
        if var_table != None:
            plot_haplotype_dissimilarity(
                legend_handles,
                y_offset,
                y_margin,
                index,
                rev_index,
                readset,
                var_table,
                genome_space,
            )

        y_offset = 0

        # Plot heatmaps
        for c_id in range(0, len(clusters)):
            if len(clusters[c_id]) < 5:
                continue
            read_id = 0
            for read in clusters[c_id]:
                start = index[readset[read][0].position]
                end = index[readset[read][-1].position]
                read_id += 1

                color_code = "C" + str(true_hap[read]) if true_hap[read] != -1 else "black"
                if color_code not in legend_handles:
                    legend_handles[color_code] = mpatches.Patch(
                        color=color_code, label=true_hap[read]
                    )
                if genome_space:
                    plt.hlines(
                        y=read_id + y_offset,
                        xmin=rev_index[start],
                        xmax=rev_index[end],
                        color=color_code,
                    )
                else:
                    plt.hlines(y=read_id + y_offset, xmin=start, xmax=end, color=color_code)

            y_offset += len(clusters[c_id]) + y_margin
            plt.hlines(y=y_offset, xmin=min_pos, xmax=max_pos, color=(0.5, 0.5, 0.5, 0.5))
            y_offset += y_margin

        plt.legend(handles=legend_handles.values(), loc="lower center", ncol=len(legend_handles))
        axes = plt.gca()
        axes.set_xlim([min_pos, max_pos])
        fig.savefig(path)
        fig.clear()

    except ImportError:
        logger.error("Plotting read clusters requires matplotlib to be installed")


def plot_haplotype_dissimilarity(
    legend_handles, y_offset, y_margin, index, rev_index, readset, var_table, genome_space=False
):
    try:
        import matplotlib

        matplotlib.use("agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        num_vars = len(readset.get_positions())
        min_pos = min(readset.get_positions()) if genome_space else 0
        max_pos = max(readset.get_positions()) if genome_space else num_vars
        readlen = avg_readlength(readset)

        # Plot heatmaps
        y_offset = 0
        y_margin = 5

        # Plot haplotype dissimilarity
        phase_vectors = get_phase(readset, var_table)
        if not phase_vectors:
            # No phasing information available
            return

        chunk = 24
        padding = int(
            readlen // 6
        )  # dissimilarity of position i is averaged over interval of 2*padding base pairs (not positions)

        if genome_space:
            # get variant density
            dens_pos = list(range(min_pos + padding, max_pos - padding, 2 * padding))
            dens_list = [
                len(list(filter(lambda x: pos - padding <= x <= pos + padding, rev_index)))
                for pos in dens_pos
            ]
            max_dens = max(dens_list)
            print("max_dens=" + str(max_dens))

            y_offset -= 104 + y_margin
            plt.hlines(y=y_offset, xmin=min_pos, xmax=max_pos, color="black", lw=1)
            plt.hlines(y=y_offset + 104, xmin=min_pos, xmax=max_pos, color="black", lw=1)
            plt.plot(
                dens_pos, [100 * x / max_dens + y_offset for x in dens_list], lw=1, color="blue"
            )

        # determines for each position, over which interval of positions the average must be taken
        intervals = []
        for i in range(num_vars):
            left = right = i
            pos = rev_index[i]
            while left - 1 >= 0 and rev_index[left - 1] >= pos - padding:
                left -= 1
            while right + 1 < num_vars and rev_index[right + 1] <= pos + padding:
                right += 1
            intervals.append([left, right])

        # One plot for each pair of haplotypes
        for i, j in it.combinations(range(len(phase_vectors)), 2):
            y_offset -= 104 + y_margin
            colors = ["C" + str(i), "C" + str(j)]
            if colors[0] not in legend_handles:
                legend_handles[colors[0]] = mpatches.Patch(color=colors[0], label=i)
            if colors[1] not in legend_handles:
                legend_handles[colors[1]] = mpatches.Patch(color=colors[1], label=j)
            dist = [
                y_offset + 2 + 100 * v
                for v in haplodist(phase_vectors[i], phase_vectors[j], intervals)
            ]
            plt.hlines(y=y_offset, xmin=min_pos, xmax=max_pos, color="black", lw=1)
            plt.hlines(y=y_offset + 104, xmin=min_pos, xmax=max_pos, color="black", lw=1)
            for k in range(ceil(num_vars / chunk)):
                start = k * chunk
                end = min(num_vars, (k + 1) * chunk + 1)
                if genome_space:
                    plt.plot(rev_index[start:end], dist[start:end], lw=1, color=colors[k % 2])
                else:
                    plt.plot(list(range(start, end)), dist[start:end], lw=1, color=colors[k % 2])

    except ImportError:
        logger.error("Plotting haplotype dissimilarities requires matplotlib to be installed")


def haplodist(h1, h2, intervals):
    if len(h1) != len(h2):
        return -1
    n = len(h1)
    return [
        relative_hamming_dist(
            h1[intervals[i][0] : min(n, intervals[i][1] + 1)],
            h2[intervals[i][0] : min(n, intervals[i][1] + 1)],
        )
        for i in range(0, n)
    ]


def relative_hamming_dist(seq1, seq2):
    if len(seq1) != len(seq2):
        return -1
    else:
        return sum([1 for i in range(len(seq1)) if seq1[i] != seq2[i]]) / len(seq1)


def draw_threading(
    readset, clustering, coverage, paths, cut_positions, haplotypes, var_table, path
):
    try:
        import matplotlib

        matplotlib.use("agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from pylab import savefig

        assert len(paths) > 0
        ploidy = len(paths[0])
        assert ploidy >= 2
        num_c = len(clustering)
        assert num_c > ploidy
        num_vars = len(coverage)

        # Detect relevant clusters
        c_map = {}
        all_threaded = set()
        for p in range(ploidy):
            for pos in range(len(paths)):
                all_threaded.add(paths[pos][p])
        c_list = sorted(all_threaded)
        for i, c_id in enumerate(c_list):
            c_map[c_id] = i
        num_c = len(c_list)

        # Setup figure
        fig = plt.figure(figsize=(num_vars / 40, num_c / 4), dpi=100)
        legend_handles = {}
        x_scale = 1
        y_margin = 0.1
        y_offset = 0
        c_height = 0.9

        # Plot cut positions
        cut_pos = cut_positions + [num_vars]
        for pos in cut_pos:
            plt.vlines(
                x=pos, ymin=0, ymax=num_c * (c_height + y_margin), color="lightgray", alpha=0.3
            )

        # Plot cluster coverage
        xs = list(range(num_vars))
        for c_id in range(num_c):
            min_pos = num_vars
            max_pos = 0
            for pos in range(num_vars):
                if c_list[c_id] in coverage[pos] and coverage[pos][c_list[c_id]] > 0:
                    min_pos = min(min_pos, pos)
                    max_pos = max(max_pos, pos)
            ys = [
                y_offset + c_height * coverage[pos][c_list[c_id]]
                if c_list[c_id] in coverage[pos]
                else y_offset
                for pos in range(min_pos, max_pos + 1)
            ]
            plt.fill_between(x=xs[min_pos : max_pos + 1], y1=ys, y2=y_offset, color="gray")
            plt.hlines(
                y=y_offset + c_height + y_margin / 2,
                xmin=0,
                xmax=x_scale * num_vars - 1,
                color="lightgray",
                alpha=0.5,
            )
            y_offset += c_height + y_margin

        # Plot paths
        for p in range(ploidy):
            legend_handles["C" + str(p)] = mpatches.Patch(
                color="C" + str(p), label=("Hap " + str(p))
            )
            current = paths[0][p]
            start = 0
            for pos in range(1, len(paths)):
                if paths[pos][p] != current:
                    plt.hlines(
                        y=(c_map[current] + 0.25 + p / ploidy * 0.5) * (c_height + y_margin),
                        xmin=x_scale * start,
                        xmax=x_scale * pos,
                        color="C" + str(p),
                        alpha=0.9,
                    )
                    current = paths[pos][p]
                    start = pos
            plt.hlines(
                y=(c_map[current] + 0.25 + p / ploidy * 0.5) * (c_height + y_margin),
                xmin=x_scale * start,
                xmax=x_scale * num_vars - 1,
                color="C" + str(p),
                alpha=0.9,
            )

        # Plot switch flip errors
        # print(cut_positions)
        # print(cut_pos)
        # print(haplotypes)

        # If we have ground truth, retrieve it
        compare = True
        try:
            phase_vectors = get_phase(readset, var_table)
            truth = []
            assert len(phase_vectors) == ploidy
            for k in range(ploidy):
                truth.append("".join(map(str, phase_vectors[k])))
        except:
            compare = False

        if compare:
            for i in range(len(cut_pos) - 1):
                block1 = [h[cut_pos[i] : min(len(paths), cut_pos[i + 1])] for h in truth]
                block2 = [h[cut_pos[i] : min(len(paths), cut_pos[i + 1])] for h in haplotypes]

                (
                    switchflips,
                    switches_in_column,
                    flips_in_column,
                    poswise_config,
                ) = compute_switch_flips_poly_bt(
                    block1,
                    block2,
                    report_error_positions=True,
                    switch_cost=1 + 1 / (num_vars * ploidy),
                )
                for pos, e in enumerate(switches_in_column):
                    if e > 0:
                        plt.vlines(
                            x=cut_pos[i] + pos,
                            ymax=-y_margin,
                            ymin=-y_margin - c_height * e,
                            color="blue",
                            alpha=0.6,
                        )
                        switches = [
                            j
                            for j in range(ploidy)
                            if poswise_config[pos][j] != poswise_config[pos - 1][j]
                        ]
                        for h in switches:
                            c_id = c_map[paths[cut_pos[i] + pos][h]]
                            plt.vlines(
                                x=cut_pos[i] + pos,
                                ymax=c_id + 0.95 * c_height,
                                ymin=c_id + 0.05 * c_height,
                                color="black",
                                alpha=0.3,
                            )
                for pos, flipped in enumerate(flips_in_column):
                    if len(flipped) == 0:
                        continue
                    if cut_pos[i] + pos >= len(paths):
                        continue
                    plt.vlines(
                        x=cut_pos[i] + pos,
                        ymax=-y_margin,
                        ymin=-y_margin - c_height * len(flipped),
                        color="orange",
                        alpha=0.6,
                    )
                    for h in flipped:
                        c_id = c_map[paths[cut_pos[i] + pos][h]]
                        plt.hlines(
                            y=(c_id + 0.25 + h / ploidy * 0.5) * (c_height + y_margin),
                            xmin=cut_pos[i] + pos - 0.5,
                            xmax=cut_pos[i] + pos + 0.5,
                            color="black",
                            alpha=0.6,
                        )

        # plt.legend(handles=legend_handles.values(), loc='lower center', ncol=len(legend_handles))
        axes = plt.gca()
        axes.set_xlim([0, num_vars - 1])
        fig.savefig(path)
        fig.clear()

    except ImportError:
        logger.error("Plotting haplotype threading requires matplotlib to be installed")
