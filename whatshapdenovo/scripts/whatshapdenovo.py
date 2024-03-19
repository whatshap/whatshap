import sys
import os
import ast
import json
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from mylog import Logger
from parallel import run_on_local, run_on_hpc
from assembly import assemble_supereads
from consensus import final_polish_parallel
from utils import *
from clustering import limitedClusterReads

__author__ = "Xiao Luo"
__version__ = "1.0.0"
__license__ = "GPL"
__date__ = 'March, 2021'

usage = """Haplotype-aware de novo assembly of diploid genome from long reads"""


def main():
    parser = ArgumentParser(prog='python whatshapdenovo.py', description=usage, epilog='Built by: {}'.
                            format(__author__), formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', dest='infile', type=str, required=True,
                        help="input file in FASTA/FASTQ format")
    parser.add_argument('-o', '--outdir', dest='outdir', type=str, required=False, default='.',
                        help="output directory")
    parser.add_argument('-t', '--threads', dest='threads', type=int, required=False, default=1,
                        help="number of threads")
    parser.add_argument('-p', '--platform', dest='platform', type=str, required=False, default='pb',
                        help="sequencing platform(PacBio CLR/PacBio HiFi/Oxford Nanopore): [pb/hifi/ont]")
    parser.add_argument('-x', '--preset', dest='preset', action='store_true', required=False,
                        help="use preset parameters")
    parser.add_argument('-g', '--genomesize', dest='genomesize', type=str, required=False, default='small',
                        help="genome size: small/large")
    parser.add_argument('--overlaps', dest='overlaps', type=str, required=False, help="input file in PAF format")
    parser.add_argument('--nsplit', dest='nsplit', type=int, required=False, default=1,
                        help="number of splitted input fasta/fastq files")
    parser.add_argument('--run_mode', dest='run_mode', type=str, required=False, default='local',
                        help="running mode: local/hpc, hpc is not available yet")
    parser.add_argument('--min_cov', dest='min_cov', type=float, required=False, default=4.0,
                        help="min coverage for trimming consensus")

    # parser.add_argument('--filter_ovlps', dest='filter_ovlps', type=ast.literal_eval, required=False,
    #                     default=True, help="filter overlaps or not, should be either True or False")

    parser.add_argument('--min_identity', dest='min_identity', type=float, required=False, default=0.75,
                        help="min identity for filtering overlaps")
    parser.add_argument('--min_read_len', dest='min_read_len', type=int, required=False, default=1000,
                        help="min read length for processing")
    parser.add_argument('--min_sread_len', dest='min_sread_len', type=int, required=False, default=1000,
                        help="min seed read length")
    parser.add_argument('--min_ovlp_len', dest='min_ovlp_len', type=int, required=False, default=1000,
                        help="min overlap length for super reads construction")

    parser.add_argument('--n_correct', dest='n_correct', type=int, required=False, default=0,
                        help="times for self error correction of raw reads")
    parser.add_argument('--n_polish', dest='n_polish', type=int, required=False, default=2,
                        help="times for super reads polishing")

    parser.add_argument('--polish_tool', dest='polish_tool', type=str, required=False, default='racon',
                        help="tool used for polishing contigs: [racon/hypo]. Note: hypo is recommended if using HiFi reads")

    parser.add_argument('--max_oh', dest='max_oh', type=int, required=False, default=1000,
                        help="max overhang length")
    parser.add_argument('--oh_ratio', dest='oh_ratio', type=float, required=False, default=0.8,
                        help="max overhang to mapping length ratio")

    parser.add_argument('--sp_oh', dest='sp_oh', type=int, required=False, default=10,
                        help="max overhang length")
    parser.add_argument('--sp_ohratio', dest='sp_ohratio', type=float, required=False, default=0.8,
                        help="max overhang to mapping length ratio")

    parser.add_argument('--sp_min_identity', dest='sp_min_identity', type=float, required=False, default=0.98,
                        help="super reads min identity for filtering overlaps")
    parser.add_argument('--sp_min_ovlplen', dest='sp_min_ovlplen', type=int, required=False, default=1000,
                        help="min overlap length for super reads construction")

    parser.add_argument('--max_tip_len', dest='max_tip_len', type=int, required=False, default=1000,
                        help="max length to be removed as tips")

    parser.add_argument('--min_cluster_size', dest='min_cluster_size', type=int, required=False, default=4,
                        help="min size of read clusters")

    parser.add_argument('--level', dest='level', type=int, required=False, default=1,
                        help="iteration times of clustering read")
    parser.add_argument('--rm_trans', dest='rm_trans', type=int, required=False, default=1,
                        help="choose to (0) keep all edges, (1) remove transitive edges, (2) to remove double transitive edges")
    parser.add_argument('--trim_ends', dest='trim_ends', type=ast.literal_eval, required=False, default=False,
                        help="trim the erroneous bases in both ends, should be either True or False")
    parser.add_argument('--rename', dest='rename', type=ast.literal_eval, required=False, default=False,
                        help="rename read name or not, should be either True or False")
    parser.add_argument('--qc', dest='qc', type=ast.literal_eval, required=False, default=False,
                        help="quality control for input reads or not, should be either True or False, TODO")

    parser.add_argument('--ctg_asm', dest='ctg_asm', type=str, required=False, default='rb',
                        help="method to assemble super reads: [rb/naive/iter], rb is recommended")

    parser.add_argument('--correct_mode', dest='correct_mode', type=str, required=False, default='msa',
                        help="method to correct raw reads: [msa/hybrid], msa is much faster than hybrid, which is recommended for large genomes")

    parser.add_argument('--max_het_snps', dest='max_het_snps', type=int, required=False, default=0,
                        help="maximum number of heterozygous SNPs to determine the contig overlap is from the identical haplotype or not")
    parser.add_argument('--min_allele_cov', dest='min_allele_cov', type=int, required=False, default=4,
                        help="number of observations of each allele")

    parser.add_argument('--n_final_polish', dest='n_final_polish', type=int, required=False, default=1,
                        help="polish times for final contigs")

    parser.add_argument('--sort_by_len', dest='sort_by_len', type=ast.literal_eval, required=False, default=True,
                        help="sort super reads by read length or not(by number of overlaps), should be either True or False")
    parser.add_argument('--rm_tmp', dest='rm_tmp', type=ast.literal_eval, required=False, default=True,
                        help="remove temp files or not, should be either True or False")
    parser.add_argument('--max_memory', dest='max_memory', type=str, required=False, default='50G',
                        help="max memory to use")
    parser.add_argument('--superead_id', dest='superead_id', type=int, required=False,
                        help="superead id required if running on HPC")

    parser.add_argument('--limited_times', dest='limited_times', type=int, required=False, default=5,
                        help="max times used for read")
    parser.add_argument('--max_ovlps', dest='max_ovlps', type=int, required=False, default=10000,
                        help="max number of overlaps for read")
    parser.add_argument('--max_cluster_size', dest='max_cluster_size', type=int, required=False, default=100000,
                        help="max cluster size")
    parser.add_argument('--version', '-v', action='version', version='%(prog)s version: 1.0.0', help='show the version')

    args = parser.parse_args()

    ## configuring log info
    log = Logger('whatshapdenovo.log', level='debug')
    # log.logger.debug('debug')
    # log.logger.info('info')
    # log.logger.warning('warning')
    default_min_cov = 4.0
    default_min_identity = 0.75
    default_n_correct = 0
    default_n_polish = 2
    default_min_cluster_size = 4
    default_level = 1
    default_sort_by_len = True
    default_trim_ends = False
    default_sp_min_ovlplen = 1000
    default_sp_min_identity = 0.98
    default_correct_mode = 'msa'
    default_ctg_asm = 'rb'
    default_n_final_polish = 1
    default_max_het_snps = 0
    default_min_allele_cov = 4

    # use preset parameters
    if args.preset:
        print('use preset parameters...')
        if args.platform == 'hifi':
            if args.genomesize == 'small':
                if args.min_cov == default_min_cov: args.min_cov = 4
                if args.min_identity == default_min_identity: args.min_identity = 0.95
                if args.n_correct == default_n_correct: args.n_correct = 0
                if args.n_polish == default_n_polish: args.n_polish = 1
                if args.min_cluster_size == default_min_cluster_size: args.min_cluster_size = 4
                if args.level == default_level: args.level = 8
                if args.sort_by_len == default_sort_by_len: args.sort_by_len = False
                if args.trim_ends == default_trim_ends: args.trim_ends = False
                if args.sp_min_identity == default_sp_min_identity: args.sp_min_identity = 0.99
                if args.correct_mode == default_correct_mode: args.correct_mode = 'hybrid'
                if args.ctg_asm == default_ctg_asm: args.ctg_asm = 'rb'
                # args.max_het_snps=0
                # args.min_allele_cov =4
                # args.n_final_polish=1
            elif args.genomesize == 'large':
                if args.min_cov == default_min_cov: args.min_cov = 4
                if args.min_identity == default_min_identity: args.min_identity = 0.0
                if args.n_correct == default_n_correct: args.n_correct = 0
                if args.n_polish == default_n_polish: args.n_polish = 1
                if args.min_cluster_size == default_min_cluster_size: args.min_cluster_size = 4
                if args.level == default_level: args.level = 8
                if args.sort_by_len == default_sort_by_len: args.sort_by_len = False
                if args.trim_ends == default_trim_ends: args.trim_ends = True
                if args.sp_min_identity == default_sp_min_identity: args.sp_min_identity = 0.98
                if args.correct_mode == default_correct_mode: args.correct_mode = 'msa'
                if args.ctg_asm == default_ctg_asm: args.ctg_asm = 'naive'
                # args.max_het_snps = 0
                # args.min_allele_cov = 4
                # args.n_final_polish = 1
            else:
                raise Exception('invalid setting for --genomesize')
        elif args.platform == 'pb':  # PacBio CLR
            if args.genomesize == 'small':
                if args.min_cov == default_min_cov: args.min_cov = 10
                if args.min_identity == default_min_identity: args.min_identity = 0.75
                if args.n_correct == default_n_correct: args.n_correct = 2
                if args.n_polish == default_n_polish: args.n_polish = 2
                if args.min_cluster_size == default_min_cluster_size: args.min_cluster_size = 10
                if args.level == default_level: args.level = 1
                if args.sort_by_len == default_sort_by_len: args.sort_by_len = True
                if args.trim_ends == default_trim_ends: args.trim_ends = True
                if args.sp_min_identity == default_sp_min_identity: args.sp_min_identity = 0.98
                if args.correct_mode == default_correct_mode: args.correct_mode = 'hybrid'
                # args.sp_oh =500
                if args.ctg_asm == default_ctg_asm: args.ctg_asm = 'rb'
                if args.max_het_snps == default_max_het_snps: args.max_het_snps = 0
                if args.min_allele_cov == default_min_allele_cov: args.min_allele_cov = 6
                if args.n_final_polish == default_n_final_polish: args.n_final_polish = 2
            elif args.genomesize == 'large':
                if args.min_cov == default_min_cov: args.min_cov = 10
                if args.min_identity == default_min_identity: args.min_identity = 0.
                if args.n_correct == default_n_correct: args.n_correct = 2
                if args.n_polish == default_n_polish: args.n_polish = 2
                if args.min_cluster_size == default_min_cluster_size: args.min_cluster_size = 10
                if args.level == default_level: args.level = 2
                if args.sort_by_len == default_sort_by_len: args.sort_by_len = True
                if args.trim_ends == default_trim_ends: args.trim_ends = True
                if args.sp_min_ovlplen == default_sp_min_ovlplen: args.sp_min_ovlplen = 5000
                if args.sp_min_identity == default_sp_min_identity: args.sp_min_identity = 0.98
                if args.correct_mode == default_correct_mode: args.correct_mode = 'msa'
                if args.ctg_asm == default_ctg_asm: args.ctg_asm = 'naive'
                if args.n_final_polish == default_n_final_polish: args.n_final_polish = 1
            else:
                raise Exception('invalid setting for --genomesize')
        elif args.platform == 'ont':
            if args.genomesize == 'small':
                if args.min_cov == default_min_cov: args.min_cov = 10
                if args.min_identity == default_min_identity: args.min_identity = 0.75
                if args.n_correct == default_n_correct: args.n_correct = 2
                if args.n_polish == default_n_polish: args.n_polish = 2
                if args.min_cluster_size == default_min_cluster_size: args.min_cluster_size = 10
                if args.level == default_level: args.level = 1
                if args.sort_by_len == default_sort_by_len: args.sort_by_len = True
                if args.trim_ends == default_trim_ends: args.trim_ends = True
                if args.sp_min_identity == default_sp_min_identity: args.sp_min_identity = 0.98
                if args.correct_mode == default_correct_mode: args.correct_mode = 'hybrid'
                # args.sp_oh = 500
                if args.ctg_asm == default_ctg_asm: args.ctg_asm = 'rb'
                if args.max_het_snps == default_max_het_snps: args.max_het_snps = 0
                if args.min_allele_cov == default_min_allele_cov: args.min_allele_cov = 6
                if args.n_final_polish == default_n_final_polish: args.n_final_polish = 2
            elif args.genomesize == 'large':
                if args.min_cov == default_min_cov: args.min_cov = 10
                if args.min_identity == default_min_identity: args.min_identity = 0.
                if args.n_correct == default_n_correct: args.n_correct = 2
                if args.n_polish == default_n_polish: args.n_polish = 2
                if args.min_cluster_size == default_min_cluster_size: args.min_cluster_size = 10
                if args.level == default_level: args.level = 2
                if args.sort_by_len == default_sort_by_len: args.sort_by_len = True
                if args.trim_ends == default_trim_ends: args.trim_ends = True
                if args.sp_min_ovlplen == default_sp_min_ovlplen: args.sp_min_ovlplen = 5000
                if args.sp_min_identity == default_sp_min_identity: args.sp_min_identity = 0.98
                if args.correct_mode == default_correct_mode: args.correct_mode = 'msa'
                if args.ctg_asm == default_ctg_asm: args.ctg_asm = 'naive'
                if args.n_final_polish == default_n_final_polish: args.n_final_polish = 1
            else:
                raise Exception('invalid setting for --genomesize')
        else:
            raise Exception('invalid setting for --platform')
    else:
        pass

    log.logger.info('splitting input fastx file into {} subfiles...'.format(args.nsplit))
    os.system("rm -rf {}/1.split_fastx".format(args.outdir))
    os.system("mkdir {}/1.split_fastx".format(args.outdir))
    # fastx_files = []
    # fastx_files.append(args.infile)
    fastx_files = preprocess_fastx(args.infile, args.outdir, args.nsplit, qc=args.qc,
                                   rename=args.rename)

    log.logger.info('splitting finished.\n')

    ovlp_files = []
    if args.overlaps:
        # use the input overlaps
        ovlp_files.append(args.overlaps)
    else:
        # compute overlaps using minimap2(no base level alignment) and filter overlaps using fpa
        log.logger.info('computing all-vs-all read overlaps...')
        os.system("rm -rf {}/2.overlap".format(args.outdir))
        os.system("mkdir {}/2.overlap".format(args.outdir))
        ovlp_files = compute_ovlps(fastx_files, args.outdir, args.threads, args.platform, args.genomesize,
                                   args.min_ovlp_len,
                                   args.min_identity, args.max_oh, args.oh_ratio)

    log.logger.info('computing overlaps finished.\n')
    log.logger.info('sorting overlaps by overlap length and matched length...')
    os.system("sort -k 11 -k 10 -nr -S {} --parallel {} {} -o {}".
              format(args.max_memory, args.threads, ovlp_files[0], ovlp_files[0]))

    log.logger.info('clustering reads...')
    # clusters_file = cluster_reads(ovlp_files, args.outdir, args.sort_by_len, args.min_cluster_size, args.level)

    clusters_file = limitedClusterReads(ovlp_files, args.outdir, args.sort_by_len, args.min_sread_len,
                                        args.min_cluster_size, args.level,
                                        args.limited_times, args.max_ovlps, args.max_cluster_size)

    # clusters_file = args.outdir + '/clustered_reads.list'

    log.logger.info('generating read and overlap files for each cluster...')
    cluster_dir = "{}/3.cluster".format(args.outdir)
    os.system('rm -rf {}'.format(cluster_dir))
    os.mkdir(cluster_dir)
    num_clusters = split_infiles_by_cluster(fastx_files, ovlp_files, clusters_file, cluster_dir, args.threads)

    # if args.rm_tmp:
    #     os.system("rm -rf {}/1.split_fastx/*".format(args.outdir))
    #     os.system("rm -rf {}/2.overlap/*".format(args.outdir))

    log.logger.info('stage1 assembly: assemble raw reads for each cluster...')
    if args.run_mode == 'local':
        run_on_local(num_clusters, args.threads, cluster_dir, args.platform, args.min_cov, args.max_tip_len,
                     args.n_correct, args.n_polish, args.rm_trans, args.trim_ends, args.polish_tool, args.rm_tmp,
                     args.correct_mode)
    elif args.run_mode == 'hpc':
        run_on_hpc()
    else:
        raise ArgumentParser("invalid running mode, must be: local or hpc", args.run_mode)

    log.logger.info('concatenate super reads sequences for all clusters...')
    all_supereads_fa = "{}/all.supereads.fa".format(args.outdir)
    os.system("for i in {}/3.cluster/c*/*.supereads.fa;do cat $i ;done >{}".
              format(args.outdir, all_supereads_fa))

    log.logger.info('stage1 assembly has been finished.')

    # do super reads assembly
    if args.ctg_asm == 'iter':  # assemble super reads iteratively
        log.logger.info('Using {} method for super reads assembly...'.format(args.ctg_asm))
        # assemble_supereads_iter()#TODO

    elif args.ctg_asm == 'rb' or args.ctg_asm == 'naive':
        log.logger.info('Using {} method for super reads assembly...'.format(args.ctg_asm))
        asm_supereads_dir = "{}/4.asm_supereads".format(args.outdir)
        os.system("rm -rf {}".format(asm_supereads_dir))
        os.system("mkdir {}".format(asm_supereads_dir))
        ctg_file, utg2supereads_old = assemble_supereads(all_supereads_fa, asm_supereads_dir, args.threads,
                                                         args.min_read_len, args.sp_min_ovlplen,
                                                         args.sp_min_identity, args.sp_oh, args.sp_ohratio,
                                                         args.max_tip_len, args.ctg_asm,
                                                         args.max_het_snps, args.min_allele_cov, args.platform,
                                                         args.rm_tmp)

        if args.n_final_polish > 0:
            log.logger.info('Polishing final contigs for {} times...'.format(args.n_final_polish))
            polish_dir = "{}/5.polish".format(args.outdir)
            os.system("rm -rf {}".format(polish_dir))
            os.system("mkdir {}".format(polish_dir))

            # rename
            utg2supereads = {}
            new2raw = {}
            with open(args.outdir + '/4.asm_supereads/reads.id_map') as fid:
                for line in fid:
                    rawid, newid = line.strip().split()
                    new2raw[newid] = rawid
            for utg, new_ids in utg2supereads_old.items():
                utg2supereads[utg] = ' '.join([new2raw[newid] for newid in new_ids])

            with open(polish_dir + '/utg2supereads.json', 'w') as fjson:
                json.dump(utg2supereads, fjson)

            final_polish_parallel(ctg_file, utg2supereads, args.polish_tool, args.n_final_polish, polish_dir,
                                  args.platform,
                                  args.threads, args.min_ovlp_len, args.min_identity)
    else:
        raise ArgumentParser('invalid method for super reads assembly, must be: rb/naive/iter', args.ctg_asm)

    os.system('mv {}/5.polish/final_contigs.polished.fa {}/contigs.fa'.format(args.outdir, args.outdir))

    if args.rm_tmp:
        os.system("rm -rf {}/3.cluster/*".format(args.outdir))
        os.system("rm -rf {}/5.polish/*".format(args.outdir))

    log.logger.info('All has been finished successfully.\n')
    log.logger.info('The final output haplotype aware contigs are here: {}/contigs.fa\n'.format(args.outdir))
    log.logger.info('Thank you for using WhatsHap-DeNovo!!\n')


if __name__ == '__main__':
    sys.exit(main())
