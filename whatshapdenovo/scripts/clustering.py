import json
import sys
import os


def limitedClusterReads(ovlp_files, outdir, sort_by_len, min_sread_len, min_cluster_size, k, limited_times, max_ovlps,
                        max_cluster_size):  # TODO: sort overlap file by overlap length,10,11 at first !!
    '''
    assign each read to limited clusters
    '''
    read2neigb = {}  # read->[read_len,neighbor reads]
    for ovlp_file in ovlp_files:
        print(ovlp_file)
        with open(ovlp_file, 'r') as fr:
            for line in fr:
                a = line.split()
                qname, qlen, tname, tlen = a[0], int(a[1]), a[5], int(a[6])
                if qname in read2neigb and len(read2neigb[qname][1]) <= max_ovlps:
                    read2neigb[qname][1].append(tname)
                else:
                    read2neigb[qname] = [qlen, [tname]]

                if tname in read2neigb and len(read2neigb[qname][1]) <= max_ovlps:
                    read2neigb[tname][1].append(qname)
                else:
                    read2neigb[tname] = [tlen, [qname]]

    print('reading overlap done...')
    # clustering through key-value
    # log.logger.info('clustering reads...')
    used_reads = {read: 0 for read in read2neigb.keys()}
    clusters_file = outdir + "/clustered_reads.list"
    fw = open(clusters_file, 'w')
    if sort_by_len:
        # TODO: Really need to sort?
        for item in sorted(read2neigb.items(), key=lambda d: d[1][0], reverse=True):
            sread = item[0]
            if used_reads[sread] or read2neigb[sread][0] < min_sread_len:
                continue
            clustered_reads = {}
            clustered_reads[sread] = 1

            last_reads = {}
            for i in range(k):
                if i == 0:
                    for r in read2neigb[sread][1]:
                        if used_reads[r] >= limited_times:
                            continue
                        last_reads[r] = 1
                        clustered_reads[r] = 1
                else:
                    tmp_reads = {}
                    if len(last_reads):
                        for sread_tmp in last_reads.keys():
                            if used_reads[sread_tmp] >= limited_times:
                                continue
                            for r in read2neigb[sread_tmp][1]:
                                if used_reads[r] >= limited_times:
                                    continue
                                clustered_reads[r] = 1
                                tmp_reads[r] = 1
                            if len(clustered_reads) > max_cluster_size:
                                break
                        last_reads = tmp_reads
                if len(clustered_reads) > max_cluster_size:
                    break
            if len(clustered_reads) < min_cluster_size:
                continue

            # mark as used reads
            for read in clustered_reads.keys():
                used_reads[read] += 1
            fw.write(' '.join([r for r in clustered_reads.keys()]) + '\n')

    else:
        for sread in read2neigb.keys():
            if used_reads[sread] or read2neigb[sread][0] < min_sread_len:
                continue
            clustered_reads = {}
            clustered_reads[sread] = 1

            last_reads = {}
            for i in range(k):
                if i == 0:
                    for r in read2neigb[sread][1]:
                        if used_reads[r] >= limited_times:
                            continue
                        last_reads[r] = 1
                        clustered_reads[r] = 1
                else:
                    tmp_reads = {}
                    if len(last_reads):
                        for sread_tmp in last_reads.keys():
                            if used_reads[sread_tmp] >= limited_times:
                                continue
                            for r in read2neigb[sread_tmp][1]:
                                if used_reads[r] >= limited_times:
                                    continue
                                clustered_reads[r] = 1
                                tmp_reads[r] = 1
                            if len(clustered_reads) > max_cluster_size:
                                break
                        last_reads = tmp_reads
                if len(clustered_reads) > max_cluster_size:
                    break
            if len(clustered_reads) < min_cluster_size:
                continue

            # mark as used reads
            for read in clustered_reads.keys():
                used_reads[read] += 1
            fw.write(' '.join([r for r in clustered_reads.keys()]) + '\n')

    del read2neigb, used_reads
    return clusters_file


if __name__ == '__main__':
    ovlp_files, outdir, sort_by_len, min_cluster_size, k = sys.argv[1:]
    max_ovlps = 90
    max_cluster_size = 30000000000000000
    limited_times = 10000000000000
    limitedClusterReads([ovlp_files], outdir, False, int(min_cluster_size), int(k), int(limited_times), max_ovlps,
                        max_cluster_size)
