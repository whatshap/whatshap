import os
import sys

def consent(id, outdir, rounds=1, type="pb"):
    '''
    correct sequencing error of raw reads using consent
    '''
    prefix = str(id)
    in_fa = outdir + '/' + prefix + ".fa"
    out_fa = ''
    logfile = "{}/{}.log".format(outdir, id)
    selfpath = sys.path[0]

    # run CONSENT-correct for x rounds
    type2 = ''
    if type == "pb" or type == "hifi":
        type2 = "PB"
    elif type == "ont":
        type2 = "ONT"
    else:
        raise Exception('invalid sequencing platform:{}'.format(type))
    tmp_fa = in_fa
    for i in range(rounds):
        out_fa = outdir + '/' + prefix + ".correct_" + str(i + 1) + ".fa"
        cmd=selfpath+"/../tools/CONSENT/CONSENT-correct -j 1 --in " + tmp_fa + " --out " + out_fa + " --type " + type2 + \
            " --tmpdir " + outdir  + " >{} 2>&1".format(logfile)# -m:200
        print('Run CONSENT correction module command:\n{}'.format(cmd))
        try:
            os.system(cmd)
        except:
            raise Exception('Error in running CONSENT, please check if CONSENT is available')
        tmp_fa = out_fa

    corrected_fa = outdir + '/' + prefix + ".corrected.fa"
    # os.system("ln -fs {} {}".format(out_fa.split('/')[-1], in_fa))
    os.system("mv {} {}".format(out_fa, corrected_fa))
    return corrected_fa


def correct_error_reads(id, outdir, rounds=1, type="pb",correct_mode='msa'):
    '''
    correct sequencing error of raw reads using consent
    '''
    corrected_fa=''
    if correct_mode=='hybrid': #msa + De Bruijn graphs
        try:
            corrected_fa=consent(id, outdir, rounds, type)
        except:
            raise Exception('Failed to perform hybrid error correction !!')

    elif correct_mode=='msa':
        prefix = str(id)
        in_fa =  prefix + ".fa"
        out_fa = ''
        logfile = "{}/{}.tmp.log".format(outdir, id)
        selfpath = sys.path[0]
        tmp_fa = in_fa
        base_dir=os.getcwd()
        os.chdir(outdir) #change work directory
        for i in range(rounds):
            if type=='pb': #Use MECAT2
                out_fa =  prefix + ".correct_" + str(i + 1) + ".fa"
                #generate config file
                config='config.mecat2.txt'
                with open(selfpath+'/../tools/MECAT2/config.template') as fr:
                    config_temp=fr.read()
                with open(config,'w') as fw:
                    fw.write("RAWREADS={}\nOUTDIR=./\n".format(tmp_fa)+config_temp)

                # correct raw reads for x rounds
                os.system('mkdir -p mecat/1-consensus/cns_pm_dir')
                os.system(selfpath+'/../tools/MECAT2/Linux-amd64/bin/mecat.pl correct {}'.format(config))
                # os.system('cat mecat/1-consensus/cns_cns_dir/p*.cns.fasta >{}'.format(out_fa)) #old MECAT2
                os.system('cat mecat/1-consensus/cns_reads.fasta >{}'.format(out_fa))
                tmp_fa=out_fa
                os.system('rm -rf mecat')

            elif type=='ont': #Use NECAT
                out_fa = prefix + ".correct_" + str(i + 1) + ".fa"
                # generate config file
                config = 'config.mecat2.txt'
                with open(selfpath + '/../tools/NECAT/config.template') as fr:
                    config_temp = fr.read()
                with open('read_list.txt','w') as fw:
                    fw.write(tmp_fa)
                with open(config, 'w') as fw:
                    fw.write("ONT_READ_LIST=read_list.txt\n" + config_temp)

                # correct raw reads for x rounds
                os.system('mkdir -p mecat/1-consensus/cns_pm_dir')
                os.system(selfpath + '/../tools/NECAT/Linux-amd64/bin/necat.pl correct {}'.format(config))
                os.system('zcat mecat/1-consensus/cns_iter1/cns.fasta.gz >{}'.format(out_fa))
                tmp_fa = out_fa
                os.system('rm -rf mecat')

            else:
                raise Exception('invalid sequencing platform:{}'.format(type))

        os.system("mv {} {}".format(out_fa, prefix + ".corrected.fa"))
        corrected_fa = outdir + '/' + prefix + ".corrected.fa"
        os.chdir(base_dir)
    return corrected_fa

def polish_seq(i, ref, reads_fa, outdir, rounds, type, polish_tool):
    '''
    polish sequence using raw or corrected reads
    :param i: the i pivot read
    :param ref: fasta which needs to be polished
    :param reads_fa: reads used to polish
    :return a polished fasta file
    '''
    prefix = str(i)
    polish_paf = outdir + '/' + prefix + '.polish.paf'
    polished_fa = outdir + '/' + prefix + '.ref.polished.fa'
    tmp_fa = outdir + '/' + prefix + '.tmp.fa'
    os.system("cp {} {}".format(ref, tmp_fa))
    logfile = "{}/{}.log".format(outdir, i)

    if polish_tool == 'racon':
        for k in range(rounds):
            if type == 'pb' or type == 'ont':
                os.system("minimap2 --secondary=no -x map-{} -c -t 1 {} {} 2>/dev/null |cut -f 1-12 >{}"
                          .format(type, tmp_fa, reads_fa, polish_paf))
            elif type == 'hifi':
                os.system("minimap2 --secondary=no -x asm20  -c -t 1 {} {} 2>/dev/null |cut -f 1-12 >{}".
                          format(tmp_fa, reads_fa, polish_paf))
            try:
                os.system("racon  -t 1 {} {} {} >{} 2>{}".format(reads_fa, polish_paf, tmp_fa, polished_fa, logfile))
            except RuntimeError as reason:
                print("Error occurs when running Racon:\n" +
                      "The reason is {}\n".format(reason) +
                      "Maybe the cluster is not homogeneous, Racon can not support this!\n" +
                      "Try compiling Racon on each machine individually if possible!",
                      file=open("{}/{}.log".format(outdir, i), 'a'))
            os.system("cp {} {}".format(polished_fa, tmp_fa))  # the input and output fasta are identical now

    elif polish_tool == 'hypo':
        polish_bam = outdir + '/' + prefix + '.polish.bam'
        polish_sort_bam = outdir + '/' + prefix + '.polish.sort.bam'
        genome_size = ''  # Approximate size of the contig, kb
        coverage = 0  # Approximate mean coverage of the reads
        flag = True
        for k in range(rounds):
            if type == 'pb' or type == 'ont':
                os.system("minimap2 --secondary=no -ax map-{} -t {} {} {} 2>/dev/null|samtools view -Sb  - > {}".
                          format(type, 1, tmp_fa, reads_fa, polish_bam))
            elif type == 'hifi':
                os.system("minimap2 --secondary=no -ax asm20 -t {} {} {} 2>/dev/null|samtools view -Sb  - > {}".
                          format(1, tmp_fa, reads_fa, polish_bam))

            os.system("samtools sort -@{} -o {} {}".format(1, polish_sort_bam, polish_bam))
            os.system("samtools index {}".format(polish_sort_bam))

            if flag:
                # only need to calculate once
                genome_size, coverage = os.popen("samtools depth " + polish_sort_bam +
                                                 " | awk '{sum+=$3} END{print NR,sum/NR}' ").read().strip().split()
                coverage = int(float(coverage)) + 1
                genome_size = str(int(float(genome_size) / 1000)) + 'k'
                flag = False

            os.system("hypo -t {} -r {} -d {} -s {} -c {} -b {} -o {} >/dev/null".
                      format(1, reads_fa, tmp_fa, genome_size, coverage, polish_sort_bam, polished_fa))
            os.system("cp {} {}".format(polished_fa, tmp_fa))
    else:
        pass

    return polished_fa


def get_posi_by_depth(lst, min_cov, w=20):
    '''
    get the start and end positions of fragments if there are low coverage regions
    in the middle of sequence, the input should be better trimmed at both ends.
    '''
    min_cov = min_cov  # * 0.8  # becuase of unstable coverage
    min_len = 20000  # TODO
    k = int(len(lst) / w) + 1
    mean_cov = 0
    start = -1
    end = -1
    flag = False
    posi_list = []
    for i in range(k):
        mean_cov = sum([int(posi_cov.strip().split()[1]) for posi_cov in lst[i * w:(i + 1) * w]]) / w
        #         print('mean_cov:{}'.format(mean_cov))
        if mean_cov < min_cov:
            if start != -1:
                if (end - start) >= min_len:
                    posi_list.append((start, end, end - start))
                start = -1
                end = -1
                flag = False
            continue
        elif mean_cov >= min_cov and not flag:
            start = int(lst[i * w].strip().split()[0])
            end = int(lst[min(len(lst), (i + 1) * w) - 1].strip().split()[0])
            flag = True
        elif flag:

            end = int(lst[min(len(lst), (i + 1) * w) - 1].strip().split()[0])
    #             print('end:{}'.format(end))
    if flag and ((end - start) >= min_len):
        posi_list.append((start, end, end - start))  # add the last one
    return posi_list


def scan_seq_by_depth(i, hap, err_consensus_fa, reads_fa, outdir, min_cov, type, trim_ends):
    '''
    map reads back to the erroneous consensus,
    only trim bases with low coverage at both ends
    '''
    out_fa = outdir + '/../{}.{}.supereads.fa'.format(i, hap)

    # compute coverage for each base, filter secondary and supplementary alignments at first
    bam = outdir + "/" + str(i) + ".polished.bam"
    if type == 'ont' or type == 'pb':
        os.system("minimap2 -ax map-" + type + " --secondary=no -t 1 " + err_consensus_fa + " " + reads_fa + \
                  "  2>/dev/null |samtools view -hS -F 2048 -|samtools sort -@ 1 - >" + bam)
    elif type == 'hifi':
        os.system("minimap2 -ax asm20 --secondary=no -t 1 " + err_consensus_fa + " " + reads_fa + \
                  "  2>/dev/null |samtools view -hS -F 2048 -|samtools sort -@ 1 - >" + bam)

    # TODO need to consider the low coverage region in the middle of sequence
    # @@ check why many errors in the middle of some super reads:
    # some reads are assigned into a wrong haplotype group or reads in the middle are belong to
    # one haplotype and reads in both ends are belong to the other, which may cause this mistake, therefore,
    # one should split the super reads into different parts and only keep the short fragments
    # (OR only keep the longest fragment) that are satisfied with min coverage requirement.
    """
    positions = os.popen("samtools depth " + bam + "|awk '$3>=" + str(min_cov) + \
                         "'|cut -f 2|sed -n '1p;$p'").read().strip()
    if positions:
        [start, end] = positions.split("\n")  # 1 based for sam
    else:
        print("{}.{}: no sequence satisfies the requirement of min coverage".format(i, hap))
        open(out_fa, 'w').close()  # new an empty file
        return
    """

    if os.path.getsize(bam) == 0:
        print("{}.{}: no read can be aligned to the super read, skipping...".format(i, hap),
              file=open("{}/../{}.log".format(outdir, i), 'a'))
        open(out_fa, 'w').close()  # new an empty file
        return
    else:
        positions = os.popen("samtools depth " + bam + "|awk '{print NR, $0}'" + "|awk '$4>=" + str(min_cov) + \
                             "'|cut -f 1 -d ' '|sed -n '1p;$p'").read().strip()
        if positions:
            [a, b] = [int(x) for x in positions.split("\n")]  # line number
            posi_cov_list = os.popen("samtools depth " + bam + "|cut -f 2,3").read().strip().split('\n')
            posi_list = get_posi_by_depth(posi_cov_list[(a - 1):b], min_cov, w=20)  # 1 based for sam
        else:
            print("{}.{}: no sequence satisfies the requirement of min coverage".format(i, hap),
                  file=open("{}/../{}.log".format(outdir, i), 'a'))
            open(out_fa, 'w').close()  # new an empty file
            return

    with open(err_consensus_fa, "r") as fr:
        seq = fr.readline()
        seq = fr.readline().strip()

    sub_k = 1
    if trim_ends:
        for start, end, _ in posi_list:
            consensus = seq[(int(start) - 1):int(end)]
            if len(posi_list) == 1:
                head = ">c_{}_{}\n".format(i, hap)
            else:
                head = ">c_{}_{}_sub{}\n".format(i, hap, sub_k)
            with open(out_fa, 'a') as fw:
                fw.write(head + consensus + '\n')
            sub_k += 1
    else:
        # only break at the misassembled positions in the middle regions
        for start, end, _ in posi_list:
            if len(posi_list) == 1:
                consensus = seq
                head = ">c_{}_{}\n".format(i, hap)
            else:
                if sub_k == 1:
                    consensus = seq[:int(end)]
                elif sub_k == len(posi_list):
                    consensus = seq[(int(start) - 1):]
                else:
                    consensus = seq[(int(start) - 1):int(end)]
                head = ">c_{}_{}_sub{}\n".format(i, hap, sub_k)
            with open(out_fa, 'a') as fw:
                fw.write(head + consensus + '\n')
            sub_k += 1

    return out_fa
