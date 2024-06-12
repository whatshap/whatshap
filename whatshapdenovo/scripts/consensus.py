from itertools import islice
import os
import json
from multiprocessing import Pool

def reverse_comp(seq):
    base2comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'n': 'n'}
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join([base2comp[base] for base in rev_seq])
    return rev_comp_seq


def get_union_seq(seq_v, seq_w, strand_v, strand_w, b1, b2, e1, e2, l1, l2, left, right):
    '''
    get the union sequence(+) of v,w and the flank length
    left=b1-b2, right=l2-e2-(l1-e1)

    from original sequences,for example:
                   b1        e1    l1
        v -------------------------> (+)
               |***|/////////|*****|
        w      <---------------------------- (-)
                   b2        e2            l2
    '''
    if strand_w == '-':
        seq_w = reverse_comp(seq_w)
        left2 = l2 - e2
        right2 = b2
    else:
        left2 = b2
        right2 = l2 - e2
    if strand_v == '-':
        left1 = left + l1 - e1
        right1 = right + b1
    else:
        left1 = left + b1
        right1 = l1 - e1 + right

    seq_left = seq_w[:left2] if left2 >= left1 else seq_v[:left1 - left2] + seq_w[:left2]
    seq_middle = seq_w[left2:l2 - right2]
    # print(len(seq_middle))

    # do not use minus index, bug will occur when index==0
    seq_right = seq_w[l2 - right2:] if right2 >= right1 else seq_w[l2 - right2:] + seq_v[len(seq_v) + right2 - right1:]
    union_seq = seq_left + seq_middle + seq_right

    # update left,right
    left = max(0, left1 - left2)
    right = max(0, right1 - right2)

    # xx=seq_v[:left]+seq_w+seq_v[-right:]
    # print("xx:{}".format(xx))
    # print("union:{}".format(union_seq))

    return union_seq, left, right


def get_err_consensus(nodes, ovlp2record, read2seq):
    """
    assume the source node(unitig) or the first node(clique) is forward(+).

    oriented reads(it does not matter if oriented or not):
                    v  ------------------------->
                        w  <----------------------------
                               ---------------------------->
                                     <------------------------------
    consensus sequence:---------------------------------------------
    v->w->... use the next node w to update the sequence of previous node v,
    the output consensus has the same error rate with input sequences.
    """
    # print("ovlp2record keys:{}".format(ovlp2record.keys()))

    sorted_edges = []
    for i in range(len(nodes) - 1):
        sorted_edges.append((nodes[i], nodes[i + 1]))

    strand_v = '+'
    strand_w = ''
    seq_v_w = ''
    left = 0
    right = 0
    l1 = b1 = e1 = l2 = b2 = e2 = 0
    begin = True

    for edge in sorted_edges:
        v, w = edge
        if strand_v == '+':
            strand_w = '+' if ovlp2record[v + ":" + w][4] == '+' else '-'
        else:
            strand_w = '-' if ovlp2record[v + ":" + w][4] == '+' else '+'
        seq_v = read2seq[v]
        seq_w = read2seq[w]

        a = ovlp2record[v + ":" + w]
        if v == a[0]:
            l1 = int(a[1])
            b1 = int(a[2])
            e1 = int(a[3])
            l2 = int(a[6])
            b2 = int(a[7])
            e2 = int(a[8])
        elif v == a[5]:
            l1 = int(a[6])
            b1 = int(a[7])
            e1 = int(a[8])
            l2 = int(a[1])
            b2 = int(a[2])
            e2 = int(a[3])

        if begin:
            seq_v_w = seq_v
            begin = False
        seq_v_w, left, right = get_union_seq(seq_v_w, seq_w, strand_v, strand_w, b1, b2, e1, e2, l1, l2, left, right)
        strand_v = strand_w

    return seq_v_w


def final_polish_xx(ctg_file, utg2supereads_tmp, polish_tool, n_final_polish, outdir, type, threads,
                 min_ovlp_len, min_identity):
    '''polish the final contigs'''

    # use the raw name for super reads, utg:0 -> [c_2_hap1,c_1_hap2,...]
    utg2supereads = {}
    new2raw = {}
    with open(outdir + '/../4.asm_supereads/reads.id_map') as fid:
        for line in fid:
            rawid, newid = line.strip().split()
            new2raw[newid] = rawid
    for utg, new_ids in utg2supereads_tmp.items():
        utg2supereads[utg] = ' '.join([new2raw[newid] for newid in new_ids])

    with open(outdir + '/utg2supereads.json', 'w') as fjson:
        json.dump(utg2supereads, fjson)

    polished_ctg_file = outdir + '/final_contigs.polished.fa'
    with open(ctg_file) as fctg:
        while True:
            ctg_item = list(islice(fctg, 2))
            if not ctg_item:
                break
            ref = outdir + '/ref.fa'
            with open(ref, 'w') as fw:
                fw.write(''.join(ctg_item))

            reads_fa = outdir + '/reads.fa'

            # write the raw reads fasta(remove duplicated reads)
            id2seq = {}
            for superead in utg2supereads[ctg_item[0].strip()[1:]].split():
                # print('supereads:{}'.format(superead))
                cluster_id, hap_id = superead.strip().split('_')[1:3]
                # it will use raw fasta if no read error correction by CONSENT
                hap_fa = outdir + '/../3.cluster/c{}/{}/{}.corrected.fa'.format(cluster_id, hap_id, cluster_id)
                with open(hap_fa) as fr:
                    while True:
                        fa_item = list(islice(fr, 2))
                        if not fa_item:
                            break
                        read_id = fa_item[0].strip()[1:]
                        read_seq = fa_item[1].strip()
                        id2seq[read_id] = read_seq
            with open(reads_fa, 'w') as fw:
                fw.write('\n'.join(['>' + id + '\n' + seq for id, seq in id2seq.items()]) + '\n')

            polished_fa = outdir + '/ref.polished.fa'
            tmp_fa = outdir + '/tmp.fa'
            os.system("cp {} {}".format(ref, tmp_fa))

            if polish_tool == 'racon':
                polish_paf = outdir + '/polish.paf'
                for k in range(n_final_polish):
                    if type == 'pb' or type == 'ont':
                        os.system(
                            "minimap2 --secondary=no -c -x map-{} -t {} {} {} 2>/dev/null|cut -f 1-12 |awk '$11>={} && $10/$11>={}' >{}".
                            format(type, threads, tmp_fa, reads_fa, min_ovlp_len, min_identity, polish_paf))
                    elif type == 'hifi':
                        os.system(
                            "minimap2 --secondary=no -c -x asm20 -t {} {} {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11>={}' >{}".
                            format(threads, tmp_fa, reads_fa, min_ovlp_len, min_identity, polish_paf))

                    os.system("racon  -t {} {} {} {} >{} 2>/dev/null".
                              format(threads, reads_fa, polish_paf, tmp_fa, polished_fa))
                    os.system("cp {} {}".format(polished_fa, tmp_fa))  # the input and output fasta are identical now
                os.system('cat {} >>{}'.format(polished_fa, polished_ctg_file))

            elif polish_tool == 'hypo':
                polish_bam = outdir + '/polish.bam'
                polish_sort_bam = outdir + '/polish.sort.bam'
                genome_size = str(int(len(ctg_item[1]) * 1. / 1000)) + 'k'  # Approximate size of the contig, kb
                coverage = 0  # Approximate mean coverage of the reads
                flag = True
                for k in range(n_final_polish):
                    if type == 'pb' or type == 'ont':
                        os.system("minimap2 --secondary=no -ax map-{} -t {} {} {} 2>/dev/null " +
                                  "| samtools view -Sb  -F 2048  - > {}".
                                  format(type, threads, tmp_fa, reads_fa, polish_bam))
                    elif type == 'hifi':
                        os.system("minimap2 --secondary=no -ax asm20 -t {} {} {} 2>/dev/null " +
                                  "| samtools view -Sb  -F 2048 - > {}".
                                  format(threads, tmp_fa, reads_fa, polish_bam))

                    os.system("samtools sort -@{} -o {} {}".format(threads, polish_sort_bam, polish_bam))
                    os.system("samtools index {}".format(polish_sort_bam))

                    if flag:
                        coverage = os.popen("samtools depth " + polish_sort_bam +
                                            " | awk '{sum+=$3} END{print sum/NR}' ").read().strip()
                        coverage = int(float(coverage)) + 1
                        flag = False

                    os.system("hypo -t {} -r {} -d {} -s {} -c {} -b {} -o {} >/dev/null".
                              format(threads, reads_fa, tmp_fa, genome_size, coverage, polish_sort_bam, polished_fa))
                    os.system("cp {} {}".format(polished_fa, tmp_fa))

                os.system('cat {} >>{}'.format(polished_fa, polished_ctg_file))

            elif polish_tool == 'arrow':
                # PacBio, Quiver [Chin et al., 2013] and its successor Arrow [Laird Smith et al., 2016]
                # conda install -c bioconda pbgcpp
                # TODO
                pass
            elif polish_tool == 'medaka':
                # ONT, NanoPolish [Loman et al., 2015] and its successor Medaka [Nanopore Technologies, 2019]
                for k in range(n_final_polish):
                    os.system("medaka_consensus -t {} -i {} -d {} -f -o {} 1>/dev/null 2>&1".
                              format(threads,reads_fa,tmp_fa,outdir+'/medaka'))

                    polished_fa=outdir+ '/medaka/consensus.fasta'
                    os.system("cp {} {}".format(polished_fa, tmp_fa))  # the input and output fasta are identical now
                os.system('cat {} >>{}'.format(polished_fa, polished_ctg_file))
            else:
                raise Exception('invalid polishing tool was found.', polish_tool)

    return polished_ctg_file

def final_polish_single(param):
    '''polish the final contigs'''
    # use the raw name for super reads, utg:0 -> [c_2_hap1,c_1_hap2,...]
    i, ctg_item, polish_tool, n_final_polish, outdir, type, threads, \
    min_ovlp_len, min_identity,utg2supereads = param

    threads=1
    outdir = outdir + '/' + i
    os.system("mkdir -p {}".format(outdir))

    ref = outdir + '/ref.fa'
    with open(ref, 'w') as fw:
        fw.write(''.join(ctg_item))

    reads_fa = outdir + '/reads.fa'

    # write the raw reads fasta(remove duplicated reads)
    id2seq = {}
    for superead in utg2supereads[ctg_item[0].strip()[1:]].split():
        # print('supereads:{}'.format(superead))
        cluster_id, hap_id = superead.strip().split('_')[1:3]
        # it will use raw fasta if no read error correction
        hap_fa = outdir + '/../../3.cluster/c{}/{}/{}.corrected.fa'.format(cluster_id, hap_id, cluster_id)
        with open(hap_fa) as fr:
            while True:
                fa_item = list(islice(fr, 2))
                if not fa_item:
                    break
                read_id = fa_item[0].strip().split()[0][1:]
                read_seq = fa_item[1].strip()
                id2seq[read_id] = read_seq
    with open(reads_fa, 'w') as fw:
        fw.write('\n'.join(['>' + id + '\n' + seq for id, seq in id2seq.items()]) + '\n')

    polished_fa = outdir + '/ref.polished.fa'
    tmp_fa = outdir + '/tmp.fa'
    os.system("cp {} {}".format(ref, tmp_fa))
    polished_seq=''
    if polish_tool == 'racon':
        polish_paf = outdir + '/polish.paf'
        for k in range(n_final_polish):
            if type == 'pb' or type == 'ont':
                os.system(
                    "minimap2 --secondary=no -c -x map-{} -t {} {} {} 2>/dev/null|cut -f 1-12 |awk '$11>={} && $10/$11>={}' >{}".
                    format(type, threads, tmp_fa, reads_fa, min_ovlp_len, min_identity, polish_paf))
            elif type == 'hifi':
                os.system(
                    "minimap2 --secondary=no -c -x asm20 -t {} {} {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11>={}' >{}".
                    format(threads, tmp_fa, reads_fa, min_ovlp_len, min_identity, polish_paf))

            os.system("racon  -t {} {} {} {} >{} 2>/dev/null".
                      format(threads, reads_fa, polish_paf, tmp_fa, polished_fa))
            os.system("cp {} {}".format(polished_fa, tmp_fa))  # the input and output fasta are identical now
        polished_seq = os.popen("cat {}".format(polished_fa)).read()


    elif polish_tool == 'hypo':
        polish_bam = outdir + '/polish.bam'
        polish_sort_bam = outdir + '/polish.sort.bam'
        genome_size = str(int(len(ctg_item[1]) * 1. / 1000)) + 'k'  # Approximate size of the contig, kb
        coverage = 0  # Approximate mean coverage of the reads
        flag = True
        for k in range(n_final_polish):
            if type == 'pb' or type == 'ont':
                os.system("minimap2 --secondary=no -ax map-{} -t {} {} {} 2>/dev/null " +
                          "| samtools view -Sb  -F 2048  - > {}".
                          format(type, threads, tmp_fa, reads_fa, polish_bam))
            elif type == 'hifi':
                os.system("minimap2 --secondary=no -ax asm20 -t {} {} {} 2>/dev/null " +
                          "| samtools view -Sb  -F 2048 - > {}".
                          format(threads, tmp_fa, reads_fa, polish_bam))

            os.system("samtools sort -@{} -o {} {}".format(threads, polish_sort_bam, polish_bam))
            os.system("samtools index {}".format(polish_sort_bam))

            if flag:
                coverage = os.popen("samtools depth " + polish_sort_bam +
                                    " | awk '{sum+=$3} END{print sum/NR}' ").read().strip()
                coverage = int(float(coverage)) + 1
                flag = False

            os.system("hypo -t {} -r {} -d {} -s {} -c {} -b {} -o {} >/dev/null".
                      format(threads, reads_fa, tmp_fa, genome_size, coverage, polish_sort_bam, polished_fa))
            os.system("cp {} {}".format(polished_fa, tmp_fa))

        polished_seq = os.popen("cat {}".format(polished_fa)).read()

    elif polish_tool == 'arrow':
        # PacBio, Quiver [Chin et al., 2013] and its successor Arrow [Laird Smith et al., 2016]
        # conda install -c bioconda pbgcpp
        # TODO
        pass
    elif polish_tool == 'medaka':
        # ONT, NanoPolish [Loman et al., 2015] and its successor Medaka [Nanopore Technologies, 2019]
        for k in range(n_final_polish):
            os.system("medaka_consensus -t {} -i {} -d {} -f -o {} 1>/dev/null 2>&1".
                      format(threads,reads_fa,tmp_fa,outdir+'/medaka'))

            polished_fa=outdir+ '/medaka/consensus.fasta'
            os.system("cp {} {}".format(polished_fa, tmp_fa))  # the input and output fasta are identical now
        polished_seq = os.popen("cat {}".format(polished_fa)).read()
    else:
        raise Exception('invalid polishing tool was found.', polish_tool)

    # os.system("rm -rf {}".format(outdir))

    return polished_seq


def final_polish_parallel(ctg_file, utg2supereads, polish_tool, n_final_polish, outdir, type, threads,
                 min_ovlp_len, min_identity):
    '''polish the final contigs'''
    # use the raw name for super reads, utg:0 -> [c_2_hap1,c_1_hap2,...]
    print('start polishing...')
    polished_ctg_file = outdir + '/final_contigs.polished.fa'

    fw = open(polished_ctg_file, 'w')
    pool = Pool(threads)
    c=0
    with open(ctg_file) as fctg:
        k = 0
        while True:
            k += threads
            if k % threads == 8:
                print('processing the {} contig...'.format(k))
            ctg_items = list(islice(fctg, 2 * threads))
            if not ctg_items:
                break
            params = []
            for i in range(len(ctg_items) // 2):
                c+=1
                param = (str(c), ctg_items[2 * i:2 * (i + 1)], polish_tool, n_final_polish, outdir, type, threads,
                        min_ovlp_len, min_identity,utg2supereads)
                params.append(param)
            polished_ctgs = pool.map(final_polish_single, params, chunksize=1)
            polished_ctgs = [ctg for ctg in polished_ctgs if ctg]

            fw.write(''.join(polished_ctgs))

    pool.close()
    pool.join()
    fw.close()

    return polished_ctg_file
