import os
import networkx as nx

from mylog import Logger
from overlap2graph import ovlp2graph
from phasing import *
from correction import *
from toolkits import *
from consensus import get_err_consensus
from supereads import filter_ovlp_based_on_reads_parallel
from filter_ovlps import rm_extra_ovlps

# log = Logger('whatshapdenovo.log', level='debug')


######################################################
################ For Raw Reads Assembly  #############
######################################################


def rename_fa(fasta, fasta2, outdir, min_len=0):
    id_map_file = outdir + '/reads.id_map'
    i = -1
    out = []
    id_map = []
    old_id = ''
    new_id = ''
    with open(fasta) as fr:
        for line in fr:
            if line.startswith('>'):
                old_id = line[1:].split()[0]
            else:
                if len(line.strip()) >= int(min_len):
                    i += 1
                    new_id = str(i)
                    out.append('>' + new_id + '\n' + line)
                    id_map.append(old_id + '\t' + new_id)
    with open(fasta2, 'w') as fw:
        fw.write(''.join(out))
    with open(id_map_file, 'w') as fw:
        fw.write('\n'.join(id_map))
    return id_map_file


def rename_paf(paf, paf2, id_map_file, min_ovlp_len=0):
    old2new = {}
    with open(id_map_file) as fr:
        for line in fr:
            a = line.strip().split()
            old2new[a[0]] = a[1]
    out = []
    with open(paf) as fr:
        for line in fr:
            a = line.strip().split()
            if (a[0] in old2new) and (a[5] in old2new) and (int(a[10]) >= min_ovlp_len):
                a[0] = old2new[a[0]]
                a[5] = old2new[a[5]]
                out.append('\t'.join(a))
    with open(paf2, 'w') as fw:
        fw.write('\n'.join(out))
    return


def get_ctg_from_simple_path(id, nodes, ovlp2record, read2seq, outdir, as_file=True):
    '''
    get supereads or contigs from a simple path[i.e., longest path/unitig path]
    '''
    # ovlp2record = read_paf("in.paf")  # OR read overlaps from paf file
    # ovlp2record = parse_paf_str(ovlp_str)
    # read2seq = get_read2seq(fasta, 'fasta')

    ctg = get_err_consensus(nodes, ovlp2record, read2seq)

    ref_file = outdir + "/" + str(id) + ".ref.fa"
    if as_file:
        with open(ref_file, 'w') as fw:
            fw.write(">" + str(id) + "\n")
            fw.write(ctg)
            return ref_file
    else:
        return ctg


def construct_digraph(digraph_file):
    edges = []
    with open(digraph_file) as fr:
        for line in fr:
            a, b = line.strip().split()
            edges.append((a, b))
    G = nx.DiGraph()
    G.add_edges_from(edges)
    # Gr = nx.dag.transitive_reduction(G) #already removed using varialquasispecies program
    return G


def longest_path(G):
    '''
    used for merging raw reads into superead
    '''
    path = nx.dag_longest_path(G)  # ['1','3','2']
    return path


def assemble_raw_reads(id, fasta, paf, outdir, max_tip_len, rm_trans):
    # rename reads id, but this is not necessary when debug finished
    fasta2 = outdir + '/' + str(id) + '.renamed.fa'
    paf2 = outdir + '/' + str(id) + '.renamed.paf'
    id_map_file = rename_fa(fasta, fasta2, outdir)
    rename_paf(paf, paf2, id_map_file)
    digraph_file = ovlp2graph(fasta2, paf2, rm_trans, threads=1, remove_inclusions='false', rm_tips='true',
                              min_read_len=0, max_tip_len=max_tip_len)
    ovlp2record = read_paf(paf2)
    read2seq = get_read2seq(fasta2, 'fasta')
    G = construct_digraph(digraph_file)
    nodes = longest_path(G)

    ref_file = get_ctg_from_simple_path(id, nodes, ovlp2record, read2seq, outdir, as_file=True)
    return ref_file


def get_superead(param):
    '''
    get the super reads for a cluster
    '''
    [i, outdir, type, min_cov, max_tip_len, n_correct, n_polish, rm_trans, trim_ends, polish_tool, rm_tmp,correct_mode] = param
    i = str(i)
    outdir = outdir + '/c' + i

    clog = Logger(outdir + '/' + i + '.log', level='debug')

    fasta = outdir + '/' + i + '.fa'
    paf = outdir + '/' + i + '.paf'
    os.system('ln -fs {} {}'.format(i + '.fa', outdir + '/' + i + '.corrected.fa'))

    if (os.path.getsize(fasta) == 0) or (os.path.getsize(paf) == 0):
        clog.logger.error("cluster {} is empty, skipping.}".format(i))
        return

    # get ad-hoc reference for calling variants
    ref = assemble_raw_reads(i, fasta, paf, outdir, max_tip_len, rm_trans)

    # polish ad-hoc reference, which can improve phasing performance
    ref = polish_seq(i, ref, fasta, outdir, rounds=1, type=type,
                     polish_tool=polish_tool)  # TODO: necessary when HiFi ??

    # get bam
    bam = get_bam(i, ref, fasta, outdir, type)

    # get vcf
    clog.logger.info("cluster:{} variant calling started...".format(i))
    vcf = call_variant(i, bam, ref, outdir, caller="longshot")

    clog.logger.info("cluster:{} reads phasing started...".format(i))
    hap2reads = phase_reads(i, vcf, bam, ref, outdir, add_unphased=True)

    read2seq = get_read2seq(fasta, mode='fasta')
    ovlp2record = read_paf(paf)

    # min_cov = 20  # min coverage of haplotype for error correction
    # haps=["hap1","hap2","unassigned"]
    for hap in hap2reads.keys():
        if len(hap2reads[hap]) < min_cov:
            clog.logger.info('cluster:{} seed read:{} only has {} reads, which is not satisfied with min coverge'. \
                             format(i, hap, len(hap2reads[hap])))
            continue
        else:
            # print(hap2reads)

            reads = hap2reads[hap]
            # print("Reads:{}".format(','.join(reads)))

            hap_outdir = outdir + '/' + hap
            os.system('mkdir -p {}'.format(hap_outdir))

            # get super read sequence for each haplotype
            hap_fasta = get_fasta(i, reads, hap_outdir, read2seq)
            # hap_fasta_raw = get_fa_raw_reads(i, reads, hap_outdir) #TODO: from other iterations

            clog.logger.info('getting {} paf file ...'.format(hap))

            hap_paf = get_paf(i, reads, hap_outdir, ovlp2record)
            clog.logger.info('getting {} paf file finished...'.format(hap_paf))
            if (not os.path.exists(hap_paf)) or (os.path.getsize(hap_paf) == 0):
                clog.logger.warning("{} does not exist or is empty, skipping super read: c_{}_{}".format(paf, i, hap))
                continue

            hap_ref = assemble_raw_reads(i, hap_fasta, hap_paf, hap_outdir, max_tip_len, rm_trans)

            # use raw reads to polish for one time, which to avoid possible bugs(large indel) caused by consent
            # when self correcting reads for multiple times.
            if n_correct > 0:
                hap_ref = polish_seq(i, hap_ref, hap_fasta, hap_outdir, rounds=1, type=type, polish_tool=polish_tool)
                corrected_fa = correct_error_reads(i, hap_outdir, rounds=n_correct, type=type,correct_mode=correct_mode)
                polished_fa = polish_seq(i, hap_ref, corrected_fa, hap_outdir, rounds=n_polish, type=type,
                                         polish_tool=polish_tool)
                trimmed_fa = scan_seq_by_depth(i, hap, polished_fa, corrected_fa, hap_outdir, min_cov, type, trim_ends)
            else:
                polished_fa = polish_seq(i, hap_ref, hap_fasta, hap_outdir, rounds=n_polish, type=type,
                                         polish_tool=polish_tool)
                # polished_fa = polish_seq(i, hap_ref, hap_fasta_raw, hap_outdir, rounds=n_polish, type=type)
                trimmed_fa = scan_seq_by_depth(i, hap, polished_fa, hap_fasta, hap_outdir, min_cov, type, trim_ends)
                os.system("ln -fs {}.fa {}/{}.corrected.fa".format(i,hap_outdir,i))
    clog.logger.info("cluster:{} super read construction finished.".format(i))
    clog.logger.info("cluster:{} super read files are:{}/{}.hap*.supereads.fa".format(i, outdir, i))

    # remove temp files
    if rm_tmp:
        for file in os.popen('ls {} |grep -v "^hap\|supereads.fa" '.format(outdir, i)).read().strip().split():
            os.system("rm -rf " + outdir + '/' + file)
        for file in os.popen(
                'ls {}/hap*/* |grep -v "{}.fa$\|{}.corrected.fa$" '.format(outdir, i, i)).read().strip().split():
            os.system("rm -rf " + file)

    return


######################################################
################ For SuperReads Assembly #############
######################################################

def cal_supereads_overlap(fasta, outdir, threads, min_ovlp_len, min_identity, o, r):
    '''
    calculate the overlaps of supereads
    '''
    paf = outdir + '/supereads.tmp.paf'
    filtered_paf = outdir + '/supereads.paf'
    # run minimap2, no cigar,step is much less computation cost.
    # TODO: still use minimap2 to calculate overlaps??

    filter_inline = 'python ' + sys.path[0] + '/filter_ovlp_inline.py {} {} {} {} '. \
        format(min_ovlp_len, min_identity, o, r)
    minimap = "minimap2 -cx ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 -c --end-bonus 100  " + \
              "-t %s %s %s |cut -f 1-12| %s >%s" % (threads, fasta, fasta, filter_inline, paf)

    os.system(minimap)
    #remove overlaps which are out of the max_ovlps limit
    # rm_extra_ovlps(paf, filtered_paf,max_ovlps=50, rm_extra_ovlps=True)
    rm_extra_ovlps(paf, filtered_paf,max_ovlps=10, rm_extra_ovlps=True) #TODO: max_ovlps = ?
    if os.path.getsize(filtered_paf)==0:
        os.system('cp {} {}/../contigs.fa'.format(fasta, outdir))
        print('No satisfied overlap found between super reads, program finished.')
        print('The final output haplotype aware contigs are here :\n{}/../contigs.fa'.format(outdir))
        sys.exit(0)
    return filtered_paf


def identify_utg_paths(G):
    '''
    identify all unitig paths, time complexity: O(V+E)
    :param G: DAG graph
    :return: [['1','2'],['4','6']], ['3','5']
    '''
    utg_paths = []
    single_nodes = {}

    s_nodes = []
    for node in G.nodes:
        if G.in_degree(node) == 0 or G.out_degree(node) >= 2 or (G.out_degree(node) == 1 and G.in_degree(node) >= 2):
            s_nodes.append(node)

    utg_path = []
    for s_node in s_nodes:
        # assign source node
        if G.out_degree(s_node) == 0:
            single_nodes[s_node] = 1  # disconnected nodes
            continue
        elif G.in_degree(s_node) > 1 and G.out_degree(s_node) > 1:
            single_nodes[s_node] = 1  # branch exists both in in&out direction
        elif G.in_degree(s_node) != 1 and G.out_degree(s_node) == 1:
            utg_path.append(s_node)  # no branch exists in out direction
        else:
            pass

        for node in G.successors(s_node):
            child_node = node
            while True:
                if G.in_degree(child_node) == 1 and G.out_degree(child_node) == 1:
                    utg_path.append(child_node)
                    child_node = list(G.successors(child_node))[0]  # update by next node
                elif G.in_degree(child_node) == 1 and G.out_degree(child_node) != 1:
                    utg_path.append(child_node)
                    break
                elif G.in_degree(child_node) > 1:
                    if G.out_degree(child_node) == 0:
                        single_nodes[child_node] = 1  # sink node with indegree > 1
                    break
            if len(utg_path) > 1:
                utg_paths.append(utg_path)
            elif len(utg_path) == 1:
                single_nodes[utg_path[0]] = 1
            utg_path = []

    return utg_paths, list(single_nodes)


def assemble_supereads(fasta, outdir, threads, min_read_len, min_ovlp_len, min_identity, o, r,
                       max_tip_len, method, max_het_snps, min_allele_cov, type, rm_tmp):
    fasta2, paf = None, None

    if method == 'rb':
        paf = cal_supereads_overlap(fasta, outdir, threads, min_ovlp_len, min_identity, o, r)
        paf = filter_ovlp_based_on_reads_parallel(paf, outdir, threads, max_het_snps, min_allele_cov, type)

        fasta2 = outdir + '/supereads.renamed.fa'
        id_map_file = rename_fa(fasta, fasta2, outdir, min_read_len)

        paf2 = outdir + '/supereads.renamed.paf'
        rename_paf(paf, paf2, id_map_file)
        paf = paf2

    elif method == 'naive':
        fasta2 = outdir + '/supereads.renamed.fa'
        rename_fa(fasta, fasta2, outdir, min_read_len)
        paf = cal_supereads_overlap(fasta2, outdir, threads, min_ovlp_len, min_identity, o, r)
    else:
        raise Exception('Warning: check the method for contig assembly.')

    digraph_file = ovlp2graph(fasta2, paf, rm_trans=1, threads=threads, remove_inclusions='true', rm_tips='true',
                              min_read_len=min_read_len, max_tip_len=max_tip_len)
    ovlp2record = read_paf(paf)
    read2seq = get_read2seq(fasta2, 'fasta')
    G = construct_digraph(digraph_file)
    utg_paths, single_nodes = identify_utg_paths(G)
    utg2supereads = {}  # save which super-reads IDs are contained in unitig. utg:0 -> ['4','1',...]
    contigs = []
    for id, path in enumerate(utg_paths):
        assert len(path) > 1
        ctg = get_ctg_from_simple_path(id, path, ovlp2record, read2seq, outdir, as_file=False)
        contigs.append(">utg:" + str(id) + "\n" + ctg + "\n")
        utg2supereads['utg:' + str(id)] = path

    for id, node in enumerate(single_nodes):
        seq = read2seq[node]
        contigs.append(">node:" + str(id) + "\n" + seq + "\n")
        utg2supereads['node:' + str(id)] = [node]

    ctg_file_rm_isolated = outdir + '/final_contigs.rm_isolated.fa'
    with open(ctg_file_rm_isolated, 'w') as fw:
        fw.write(''.join(contigs))

    # add supreads which are not included in digraph
    # reads_in_G={ v:1 for v in G.nodes}
    id = 0
    for read in read2seq.keys():
        if read not in G:
            id += 1
            seq = read2seq[read]
            contigs.append(">isolated_node:" + str(id) + "\n" + seq + "\n")
            utg2supereads['isolated_node:' + str(id)] = [read]

    ctg_file = outdir + '/final_contigs.fa'
    with open(ctg_file, 'w') as fw:
        fw.write(''.join(contigs))

    if rm_tmp:
        os.system("rm -rf {}/supereads.*  {}/overlaps.*".format(outdir, outdir))
    return ctg_file, utg2supereads


if __name__ == '__main__':
    fasta, outdir, threads, min_read_len, min_ovlp_len, min_identity, o, r, max_tip_len = sys.argv[1:]
    assemble_supereads(fasta, outdir, threads, min_read_len, min_ovlp_len, min_identity, o, r, max_tip_len)
