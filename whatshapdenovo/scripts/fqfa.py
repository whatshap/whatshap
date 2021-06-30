import sys

'''
 Convert a fastq file to a fasta file and vice versa.
'''

def fq2fa(fq,fa):
    out=[]
    with open(fq,'r') as fr:
        c=0
        for line in fr:
            c+=1
            if (c%4) == 1:
                out.append('>' + line[1:])
            elif (c%4) == 2:
                out.append(line)
    with open(fa,'w') as fw:
        fw.write(''.join(out))
    return

def fa2fq(fa,fq):
    seq=''
    out=[]
    n_reads=0
    with open(fa,'r') as fr:
        for line in fr:
            if line.startswith('>'):
                n_reads+=1
                if seq:
                    score='I'*len(seq)
                    out.append(seq+'\n+\n'+score+'\n')
                    seq=''
                out.append('@'+line[1:])
            else:
                seq+=line.strip()
    #the last one
    score='I'*len(seq)
    out.append(seq+'\n+\n'+score+'\n')
    with open(fq,'w') as fw:
        fw.write(''.join(out))
    return n_reads

def main():

    if len(sys.argv) < 3:
        print('ERROR: please check input parameters.')
        return 1
          
    infile = sys.argv[1]
    outfile = sys.argv[2]
    mode=sys.argv[3]
    if mode =='fq2fa':
        fq2fa(infile,outfile)
    elif mode== 'fa2fq':
        fa2fq(infile,outfile)
    else:
        print('mode error')
        return 1


if __name__ == '__main__':
        sys.exit(main())
