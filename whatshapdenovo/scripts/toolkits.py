import networkx as nx
import sys

def fq_or_fa(file):
    with open(file) as fr:
        s = fr.readline()[0]
    mode = ''
    if s == '>':
        mode = 'fasta'
    elif s == '@':
        mode = 'fastq'
    else:
        raise Exception("invalid input file, must be FASTA/FASTQ format.", file)
    return mode


def get_read2seq(file, mode='fastq'):
    read2seq = {}
    read = ''
    seq = ''
    if mode == 'fastq':
        with open(file) as fr:
            for i, line in enumerate(fr):
                if i % 4 == 0:
                    read = line.strip().strip('@')
                elif i % 4 == 1:
                    seq = line.strip()
                    read2seq[read] = seq
                else:
                    continue
    elif mode == 'fasta':
        with open(file) as fr:
            for i, line in enumerate(fr):
                if i % 2 == 0:
                    read = line.strip().strip('>')
                elif i % 2 == 1:
                    seq = line.strip()
                    read2seq[read] = seq
                else:
                    continue
    else:
        print("Error: unknown mode, only fastq or fasta permitted.")
        sys.exit(1)
    return read2seq

def get_fasta(i, nodes, outdir,read2seq):
    '''
    :param nodes: reads ID list
    :return: fasta file
    '''
    prefix = str(i)
    fa = outdir + '/' + prefix + '.fa'
    with open(fa, 'w') as fw:
        for node in nodes:
            fw.write(">" + node + "\n" + read2seq[node] + "\n")
    return fa



def get_paf(i, reads, outdir,ovlp2record):
    '''
    generate paf file from given read ids
    '''
    # print(reads)
    paf = outdir + '/' + str(i) + '.paf'
    out = []
    for i in range(len(reads) - 1):
        for j in range(i + 1, len(reads)):
            key = reads[i] + ':' + reads[j]
            if key in ovlp2record:
                out.append('\t'.join(ovlp2record[key]))
            # else:
            # print('{} does not exist in overlap file!'.format(key))
    # print(list(ovlp2record.keys()))
    # print('paf out:')
    # print('\n'.join(out))
    with open(paf, 'w') as fw:
        fw.write('\n'.join(out))
    return paf


def read_paf(paf):
    ovlp2record = {}  # store the record of v2w and w2v
    with open(paf, 'r') as fr:
        for line in fr:
            a = line.strip().split()
            node1 = a[0]
            node2 = a[5]
            key = node1 + ":" + node2
            key2 = node2 + ":" + node1
            ovlp2record[key] = a
            ovlp2record[key2] = ovlp2record[key]
    return ovlp2record
