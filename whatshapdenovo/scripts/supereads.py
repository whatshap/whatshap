import os
import sys
from multiprocessing import Pool
from itertools import islice
from phasing import get_bam, call_variant
import assembly


def merge_seq4pair_reads(ctg_fa1, ctg_fa2, ovlp_item, tmpdir):
    ovlp2record = {}
    a = ovlp_item.strip().split()
    node1 = a[0]
    node2 = a[5]
    key = node1 + ":" + node2
    key2 = node2 + ":" + node1
    ovlp2record[key] = a
    ovlp2record[key2] = ovlp2record[key]

    read2seq1 = assembly.get_read2seq(ctg_fa1, 'fasta')
    read2seq2 = assembly.get_read2seq(ctg_fa2, 'fasta')
    read2seq = dict(read2seq1, **read2seq2)  # merge

    b1 = int(a[2])
    b2 = int(a[7]) if a[4] == '+' else int(a[6]) - int(a[8])
    reads = [node1, node2] if b1 > b2 else [node2, node1]
    to_direction = 'first_to_second' if b1 > b2 else 'second_to_first'
    ref = assembly.get_ctg_from_simple_path('out', reads, ovlp2record, read2seq, tmpdir, as_file=True)
    return ref, to_direction


def filter_ovlp_based_on_reads(param):
    ovlp_item, max_het_snps, min_allele_cov, outdir, type = param
    basedir = outdir + '/..'
    a = ovlp_item.split()
    c1, c1_hap = a[0].strip().split('_')[1:3]
    c2, c2_hap = a[5].strip().split('_')[1:3]

    #filter overlaps from the same cluster, like this: c_12_hap1,c_12_hap2.
    #TODO:what about c_12_hap1_sub1,c_12_hap2_sub2
    if c1 == c2:
        return
    tmpdir = outdir + '/' + a[0] + '.' + a[5]
    if os.path.exists(tmpdir):
        # this is a duplicated overlap
        return
    else:
        os.system("mkdir {}".format(tmpdir))

    raw_fa1 = basedir + '/3.cluster/c' + c1 + '/' + c1_hap + '/' + c1 + '.fa'  # hap-level cluster 1
    raw_fa2 = basedir + '/3.cluster/c' + c2 + '/' + c2_hap + '/' + c2 + '.fa'  # hap-level cluster 2
    ctg_fa1 = basedir + '/3.cluster/c' + c1 + '/' + c1 + '.' + c1_hap + '.supereads.fa'
    ctg_fa2 = basedir + '/3.cluster/c' + c2 + '/' + c2 + '.' + c2_hap + '.supereads.fa'

    # get ad-hoc ref
    i = 'out'
    ref, to_direction = merge_seq4pair_reads(ctg_fa1, ctg_fa2, ovlp_item, tmpdir)

    raw_fa=tmpdir+'/raw.fa' # merge fa1, fa2 to unique sequence
    seq_dict={}
    for fasta in [raw_fa1,raw_fa2]:
        with open(fasta) as fr:
            while True:
                lines=''.join(list(islice(fr,2)))
                if not lines:
                    break
                seq_dict[lines]=1
    with open(raw_fa,'w') as fw:
        fw.write(''.join(seq_dict.keys()))
    # get bam
    bam = get_bam(i, ref, raw_fa, tmpdir, type)

    # get vcf
    vcf = call_variant(i, bam, ref, tmpdir, caller="longshot")
    n_het_snps = 0
    if os.path.getsize(vcf):
        with open(vcf) as fr:
            for line in fr:
                if line.startswith('#'):
                    continue
                else:
                    if not line:
                        continue
                    elif line.strip().split()[-1].split(':')[0] == '0/1':
                        vcf_items=line.split()
                        # if vcf_items[6] != 'PASS':
                        #     continue
                        pos = int(vcf_items[1])
                        if to_direction == 'first_to_second':
                            ovlp_minpos = int(a[2])
                            ovlp_maxpos = int(a[3])
                        else:
                            if a[4]=='+':
                                ovlp_minpos = int(a[7])
                                ovlp_maxpos = int(a[8])
                            else:
                                ovlp_minpos = max(int(a[1])-int(a[3]),int(a[7]))
                                ovlp_maxpos = ovlp_minpos+int(a[8])-int(a[7])

                        # only consider the SNVs in the overlap region
                        print('minpos:{} \nmaxpos:{}.'.format(ovlp_minpos,ovlp_maxpos), file=open(tmpdir + '/snp.log', 'a'))
                        if pos > ovlp_maxpos or pos < ovlp_minpos:
                            continue

                        # Number of Observations of Each Allele
                        # DP=20;AC=13,6;AM=1;MC=0
                        alleles_cov =[int(i) for i in vcf_items[7].split(';')[1].replace('AC=', '').split(',')]
                        if alleles_cov[0] >= min_allele_cov and alleles_cov[1] >= min_allele_cov:
                            if max(alleles_cov)/min(alleles_cov) < 2:#TODO !!
                                n_het_snps += 1

    print('{} heterozygous SNPs are found.'.format(n_het_snps),file=open(tmpdir+'/snp.log','a'))
    os.system("rm -rf {}".format(tmpdir))

    if n_het_snps <= max_het_snps:
        return ovlp_item
    else:
        return


def filter_ovlp_based_on_reads_parallel(paf, outdir, threads, max_het_snps, min_allele_cov, type):
    '''
    filter overlaps which are from different haplotypes based on raw reads,
    filter the overlap if there is no enough evidence: the number of heterozygous variants >= threshold
    '''
    filtered_paf = outdir + '/supereads.filtered.paf'
    fw = open(filtered_paf, 'w')
    pool = Pool(threads)

    with open(paf) as fr:
        while True:
            next_n_lines = list(islice(fr, threads))
            if not next_n_lines:
                break
            params = [(line, max_het_snps, min_allele_cov, outdir, type) for line in next_n_lines]
            ovlp_list = pool.map(filter_ovlp_based_on_reads, params, chunksize=1)
            ovlp_list = [ovlp for ovlp in ovlp_list if ovlp]
            fw.write(''.join(ovlp_list))
    pool.close()
    pool.join()
    fw.close()
    if os.path.getsize(filtered_paf)==0:
        os.system('cp {}/../all.supereads.fa {}/../contigs.fa'.format(outdir, outdir))
        print('No satisfied overlap found between super reads, program finished.')
        print('The final output haplotype aware contigs are here :\n{}/../contigs.fa'.format(outdir))
        sys.exit(0)
    return filtered_paf


## deprecated functions, maybe
def get_hapsread2rawreads(outdir):
    '''
    get a dict to store the super-read of each haplotype: raw reads
    from the previous iteration output
    :param pre_outdir:
    :return:
    '''
    # outdir=prefix + "/iter"+str(n)
    fastas = os.popen("ls {}/p_*/hap*/[0-9]*.raw_reads.fa".format(outdir)).read().strip().split('\n')
    hapsread2rawreads = {}

    for fasta in fastas:
        # print("fasta:{}".format(fasta))
        a = fasta.split('/')
        hapsread = a[-3] + ":" + a[-2]  # haplotype super read
        raw_reads = []
        with open(fasta, 'r') as fr:
            for line in fr:
                if line.startswith('>'):
                    raw_reads.append(line.strip().replace('>', ''))
        hapsread2rawreads[hapsread] = raw_reads

    return hapsread2rawreads


if __name__ == '__main__':
    pass
