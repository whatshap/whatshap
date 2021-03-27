from itertools import combinations, islice
from multiprocessing import Pool
import sys
import os
import copy

from toolkits import fq_or_fa
from mylog import Logger


# log = Logger('whatshapdenovo.log', level='debug')


def process_read(param):
    next_n_lines, file_id, iter_num, r, outfile, qc, rename, mode = param
    out_lines = []
    for i, line in enumerate(next_n_lines):
        if mode == 'fasta':
            if rename:
                if i % 2 == 0:
                    read_id = str(r * (iter_num - 1) + (i // 2 + 1))
                    out_lines.append('>{}_{}\n'.format(file_id, read_id))
                elif i % 2 == 1:
                    out_lines.append(line)
            else:
                out_lines.append(line)
        elif mode == 'fastq':  # TODO: if QC
            if rename:
                if i % 4 == 0:
                    read_id = str(r * (iter_num - 1) + (i // 4 + 1))
                    out_lines.append('>{}_{}\n'.format(file_id, read_id))
                elif i % 4 == 1:
                    out_lines.append(line)
                else:
                    pass
            else:
                if i % 4 == 0:
                    out_lines.append(line.replace('@', '>'))
                elif i % 4 == 1:
                    out_lines.append(line)
                else:
                    pass
        else:
            raise Exception("invalid input file, must be FASTA/FASTQ format.")

    with open(outfile, 'a') as fw:
        fw.write(''.join(out_lines))
    return


def preprocess_fastx(infile, outdir, nsplit, qc, rename):
    fastx_files = []
    for file_id in range(nsplit):
        fastx_file = outdir + "/1.split_fastx/s" + str(file_id + 1) + ".fa"
        fastx_files.append(fastx_file)

    mode = fq_or_fa(infile)
    r = 100  # number of reads processed by single thread once a time
    n_lines = r * 4 if mode == 'fastq' else r * 2
    pool = Pool(nsplit)
    to_end = False
    iter_num = 1
    with open(infile) as fr:
        while True:
            params = []
            for i in range(nsplit):
                next_n_lines = list(islice(fr, n_lines))  # much faster than fr.readlines( n_lines)
                if not next_n_lines:
                    to_end = True
                    break
                param = (next_n_lines, str(i + 1), iter_num, r, fastx_files[i], qc, rename, mode)
                params.append(param)
            pool.map(process_read, params, chunksize=1)  # jobs are sorted
            iter_num += 1
            if to_end: break

    pool.close()
    pool.join()

    return fastx_files


def compute_ovlp(fastx_i, fastx_j, outdir, threads, platform,genomesize, min_ovlp_len, min_identity,
                 max_oh, oh_ratio):
    '''compute overlaps from a pair of fasta/fastq file
    & filter overlaps using fpa(cannot set overhang cutoff and cannot remove duplicated overlaps)
    OR my own script?
    '''
    prefix_i = fastx_i.split('/')[-1].replace('.fa', '')
    prefix_j = fastx_j.split('/')[-1].replace('.fa', '')
    ovlp_file = "{}/2.overlap/{}_{}.paf".format(outdir, prefix_i, prefix_j)
    if genomesize=='small':
        if platform == 'pb' :
            os.system("minimap2 -cx ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 -t {} \
                        {}  {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11 >={} ' |fpa drop -i -m  >{}"
                      .format(threads, fastx_i, fastx_j, min_ovlp_len, min_identity, ovlp_file)
                      )
        elif platform == 'hifi': # for HiFi reads, it seems no need to use -c, which largely speeds up.
            os.system("minimap2 -cx ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 -t {} \
                        {}  {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11 >={} ' |fpa drop -i -m  >{}"
                      .format(threads, fastx_i, fastx_j, min_ovlp_len, min_identity, ovlp_file)
                      )
        elif platform == 'ont':# use cut -f 1-12 and then use fpa to prevent big RAM.
            print('ONT platform')
            os.system("minimap2 -cx ava-ont -k15 -Xw5 -m100 -g10000 -r2000 --max-chain-skip 25  -t {} \
                        {}  {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11 >={} ' |fpa drop -i -m  >{}"
                      .format(threads, fastx_i, fastx_j, min_ovlp_len, min_identity, ovlp_file)
                      )
        else:
            raise Exception('Unvalid sequencing platform, please set one of them: pb/hifi/ont')
    else:
        if platform == 'pb' :
            os.system("minimap2 -x ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 -t {} \
                        {}  {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11 >={} ' |fpa drop -i -m  >{}"
                      .format(threads, fastx_i, fastx_j, min_ovlp_len, min_identity, ovlp_file)
                      )
        elif platform == 'hifi': # for HiFi reads, it seems no need to use -c, which largely speeds up.
            os.system("minimap2 -x ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 -t {} \
                        {}  {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11 >={} ' |fpa drop -i -m  >{}"
                      .format(threads, fastx_i, fastx_j, min_ovlp_len, min_identity, ovlp_file)
                      )
        elif platform == 'ont':
            print('ONT platform')
            os.system("minimap2 -x ava-ont -k15 -Xw5 -m100 -g10000 -r2000 --max-chain-skip 25  -t {} \
                        {}  {} 2>/dev/null |cut -f 1-12 |awk '$11>={} && $10/$11 >={} ' |fpa drop -i -m  >{}"
                      .format(threads, fastx_i, fastx_j, min_ovlp_len, min_identity, ovlp_file)
                      )
        else:
            raise Exception('Unvalid sequencing platform, please set one of them: pb/hifi/ont')
    return ovlp_file


def compute_ovlps(fastx_files, outdir, threads, platform,genomesize, min_ovlp_len, min_identity,
                  max_oh, oh_ratio):
    '''
    compute overlaps from splitted fasta/fastq file list

    '''
    ovlp_files = []
    # TODO: Minimap2 gives different result for different orders of input files.
    # and different with combined input results, so need to add x1.fa_x2.fa & x2.fa_x1.fa
    # it seems splitting input files leads less overlaps
    for fastx_i, fastx_j in list(combinations(fastx_files, 2)) + [(f, f) for f in fastx_files]:
        ovlp_file = compute_ovlp(fastx_i, fastx_j, outdir, threads, platform, genomesize,min_ovlp_len, min_identity,
                                 max_oh, oh_ratio)
        ovlp_files.append(ovlp_file)
    return ovlp_files


def cluster_reads(ovlp_files, outdir, sort_by_len, min_cluster_size, k):
    '''
    :param sort_by_len:
    :param k: iteration times of clustering through overlapped reads
    :return:
    '''
    read2neigb = {}  # read->[read_len,neighbor reads]
    for ovlp_file in ovlp_files:
        with open(ovlp_file, 'r') as fr:
            for line in fr:
                a = line.split()
                qname, qlen, tname, tlen = a[0], int(a[1]), a[5], int(a[6])
                if qname in read2neigb:
                    read2neigb[qname][1] = read2neigb[qname][1] + " " + tname
                else:
                    read2neigb[qname] = [qlen, tname]

                if tname in read2neigb:
                    read2neigb[tname][1] = read2neigb[tname][1] + " " + qname
                else:
                    read2neigb[tname] = [tlen, qname]

    # clustering through key-value
    # log.logger.info('clustering reads...')
    clusters_file = outdir + "/clustered_reads.list"
    fw = open(clusters_file, 'w')
    if sort_by_len:
        # TODO: maybe there is bug ? cannot sort since some length was set as 0. Really need to sort?
        for item in sorted(read2neigb.items(), key=lambda d: d[1][0], reverse=True):
            sread = item[0]
            if not read2neigb[sread][0]:
                continue
            clustered_reads = {}
            clustered_reads[sread] = 1

            last_reads = {}
            for i in range(k):
                if i == 0:
                    for r in read2neigb[sread][1].split():
                        last_reads[r] = 1
                        clustered_reads[r] = 1
                else:
                    tmp_reads = {}
                    for sread_tmp in last_reads.keys():
                        for r in read2neigb[sread_tmp][1].split():
                            clustered_reads[r] = 1
                            tmp_reads[r] = 1
                    last_reads = tmp_reads
            if len(clustered_reads) < min_cluster_size:
                continue
            for read in clustered_reads.keys():
                read2neigb[read][0] = 0  # mark as used reads
            fw.write(' '.join([r for r in clustered_reads.keys()]) + '\n')
    else:
        for sread in read2neigb.keys():
            if not read2neigb[sread][0]:
                continue
            clustered_reads = {}
            clustered_reads[sread] = 1

            last_reads = {}
            for i in range(k):
                if i == 0:
                    for r in read2neigb[sread][1].split():
                        last_reads[r] = 1
                        clustered_reads[r] = 1
                else:
                    tmp_reads = {}
                    for sread_tmp in last_reads.keys():
                        for r in read2neigb[sread_tmp][1].split():
                            clustered_reads[r] = 1
                            tmp_reads[r] = 1
                    last_reads = tmp_reads

            if len(clustered_reads) < min_cluster_size:
                continue
            for read in clustered_reads.keys():
                read2neigb[read][0] = 0  # mark as used reads
            fw.write(' '.join([r for r in clustered_reads.keys()]) + '\n')

    return clusters_file


def split_infiles_by_cluster(fastx_files, ovlp_files, clusters_file, outdir, threads):
    read2clusters = {}
    with open(clusters_file) as fr:
        cluster_id = 0
        for line in fr:
            cluster_id += 1
            os.system("mkdir -p {}/c{}".format(outdir, cluster_id))
            for read in line.strip().split():
                if read in read2clusters:
                    read2clusters[read] = read2clusters[read] + ' ' + str(cluster_id)
                else:
                    read2clusters[read] = str(cluster_id)
    num_clusters = cluster_id

    # generate fasta for each cluster
    cutoff = 1000000
    num_reads = 0
    cluster2fa = {}

    for fastx_file in fastx_files:
        with open(fastx_file) as fr:
            while True:
                read_item = list(islice(fr, 2))
                if not read_item:
                    break
                num_reads += 1
                read = read_item[0][1:].strip()
                if read in read2clusters:
                    for cluster_id in set(read2clusters[read].split()):
                        if cluster_id in cluster2fa:
                            cluster2fa[cluster_id].append(''.join(read_item))
                        else:
                            cluster2fa[cluster_id] = [''.join(read_item)]

                if num_reads > cutoff:
                    for cluster_id in cluster2fa.keys():
                        out_fa = outdir + "/c" + cluster_id + "/" + cluster_id + ".fa"
                        with open(out_fa, 'a') as fw:
                            fw.write(''.join(cluster2fa[cluster_id]))
                    num_reads = 0
                    cluster2fa = {}
    # write the last part
    for cluster_id in cluster2fa.keys():
        out_fa = outdir + "/c" + cluster_id + "/" + cluster_id + ".fa"
        with open(out_fa, 'a') as fw:
            fw.write(''.join(cluster2fa[cluster_id]))
    cluster2fa = {}

    # generate paf for each cluster
    cutoff = 50000000
    num_ovlps = 0
    k = 1000
    cluster2ovlps = {}
    for ovlp_file in ovlp_files:
        with open(ovlp_file) as fr:
            while True:
                next_n_lines = list(islice(fr, k))
                if not next_n_lines:
                    break
                num_ovlps += k
                for line in next_n_lines:
                    a = line.split()
                    qname, tname = a[0], a[5]
                    if (qname in read2clusters) and (tname in read2clusters):
                        for cluster_id in set(read2clusters[qname].split()).intersection(
                                set(read2clusters[tname].split())):
                            if cluster_id in cluster2ovlps:
                                cluster2ovlps[cluster_id].append(line)
                            else:
                                cluster2ovlps[cluster_id] = [line]
                if num_ovlps > cutoff:
                    for cluster_id in cluster2ovlps.keys():
                        out_paf = outdir + "/c" + cluster_id + "/" + cluster_id + ".paf"
                        with open(out_paf, 'a') as fw:
                            fw.write(''.join(cluster2ovlps[cluster_id]))
                    num_ovlps = 0
                    cluster2ovlps = {}
    # write the last part
    for cluster_id in cluster2ovlps.keys():
        out_paf = outdir + "/c" + cluster_id + "/" + cluster_id + ".paf"
        with open(out_paf, 'a') as fw:
            fw.write(''.join(cluster2ovlps[cluster_id]))

    cluster2ovlps = {}

    return num_clusters


## testing functions ##
if __name__ == '__main__':
    # preprocess_fastx(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])

    cluster_reads([sys.argv[1]], sys.argv[2], bool(sys.argv[3]), int(sys.argv[4]))
    # cluster_reads(ovlp_files, outdir, sort_by_len, k):
