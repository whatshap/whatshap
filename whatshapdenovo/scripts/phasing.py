import os
import re


def get_bam(i, ref, fasta, outdir, type):
    prefix = str(i)
    bam = outdir + '/' + prefix + ".bam"
    # TODO: is cigar necessary?
    if type == 'pb' or type == 'ont':
        os.system("minimap2 -ax map-" + type + " -t 1 --secondary=no " + ref + " " + fasta + \
                  " 2>/dev/null |samtools view -hS -F 2048 - |samtools sort -@ 1 - >" + bam)
    elif type == 'hifi':
        os.system("minimap2 -ax asm20 -t 1 --secondary=no " + ref + " " + fasta + \
                  " 2>/dev/null|samtools view -hS -F 2048 - |samtools sort -@ 1 - >" + bam)
    return bam


def call_variant(i, bam, ref, outdir, caller="longshot"):
    '''
    iteration1:call variant using longshot which is applicable only for long reads
    subsequent iterations: use bcftools(--skip-indels because of low speed)
    because longshot cannot work on low coverage super-reads
    '''
    vcf = "{}/{}.nophase.vcf".format(outdir, i)
    os.system("touch {}".format(vcf))
    logfile = "{}/{}.log".format(outdir, i)
    if caller == "longshot":
        # build index
        os.system("samtools index {}".format(bam))
        os.system("samtools faidx {}".format(ref))
        longshot_cmd = "longshot --bam {} --ref {} -n -F --out {} >{} 2>&1". \
            format(bam, ref, vcf, logfile)
        # print(longshot_cmd)
        os.system(longshot_cmd)
    elif caller == "bcftools":
        os.system("samtools index {}".format(bam))
        os.system("samtools faidx {}".format(ref))
        cmd = "bcftools mpileup --skip-indels -f {} {} |bcftools call -mv -Ov -o {}". \
            format(ref, bam, vcf)
        os.system(cmd)
    else:
        raise Exception("Invalid caller, should be longshot or bcftools.")
    return vcf


def phase_reads(i, vcf, bam, ref, outdir, add_unphased=True):
    '''
    phase reads using whatshap which is much faster than longshot(built in Hapcut2)
    at the sacrifice of accuracy.
    '''
    hap2reads = {}
    log_file = open("{}/{}.log".format(outdir, i), "a")
    # if no variant in vcf
    if os.path.getsize(vcf) == 0:
        print("cannot phase reads in this super read: c{} because of no variant".format(i),
              file=log_file)
        reads = os.popen("samtools view {}|cut -f 1".format(bam)).read().strip().split('\n')
        hap2reads['hapUnknown'] = reads
        return hap2reads

    phased_vcf = outdir + '/{}.phased.vcf'.format(i)
    if os.system("whatshap phase -o {} --reference {} --ignore-read-groups {} {} >/dev/null 2>&1".format(phased_vcf, ref, vcf, bam)):
        print('ERROR: whatshap failed...', file=log_file)
        return
    os.system("bgzip -f {} >/dev/null 2>&1".format(phased_vcf))
    os.system("tabix -p vcf -f {}.gz >/dev/null 2>&1".format(phased_vcf))
    phase_cmd = "whatshap haplotag --reference {} --ignore-read-groups {}.gz {} 2>/dev/null| samtools view -|cut -f 1,12- ". \
        format(ref, phased_vcf, bam)
    phased_info = os.popen(phase_cmd).read().strip()
    # print('phase cmd:{}'.format(phase_cmd))
    # print("phased_info:{}".format(phased_info))

    for line in phased_info.split('\n'):
        read = line.strip().split()[0]
        hap = ''
        # hap_info HP:i:1
        if line.find('HP:i:') + 1:
            d = re.match(r'.+HP:i:(\d+)\s+', line).group(1)
            hap = 'hap' + d
        else:
            hap = 'hapUnknown'

        if hap in hap2reads:
            hap2reads[hap].append(read)
        else:
            hap2reads[hap] = [read]

    # add unphased reads into hap1 & hap2 groups and remove unknown group
    if add_unphased and ('hap1' in hap2reads) and ('hap2' in hap2reads) and ('hapUnknown' in hap2reads):
        hap2reads['hap1'].extend(hap2reads['hapUnknown'])
        hap2reads['hap2'].extend(hap2reads['hapUnknown'])
        del hap2reads['hapUnknown']

    return hap2reads
