#!/usr/bin/env python

import os
import sys
import subprocess

from fqfa import fa2fq

"""
This program is used to convert overlap file (PAF format) to directed acyclic graph(DAG) file.
It is based on c++ code of Viralquasispecies executable program.

"""


def ovlp2graph(fasta, paf, rm_trans=0, threads=1, remove_inclusions='false', rm_tips='true', min_read_len=1000,
               max_tip_len=1000):
    # id in fasta and paf files must be numerical
    # Note: the read ID in fasta must be from 0 to N,and so are IDs in overlap file.
    # otherwise, ID will be converted automatically.
    # It seems there is a bug when removing single transitive edges(no bug if >1) in viralquasispecies,
    # so only remove cycles and tips.
    outdir = os.path.dirname(os.path.abspath(fasta))
    selfpath = sys.path[0]
    # min_read_len = 1000
    # max_tip_len = 10000  # default: read length
    viralquasispecies = selfpath + "/../bin/ViralQuasispecies"
    verbose = 'true'

    single_fq = outdir + '/singles.fastq'
    n_reads = fa2fq(fasta, single_fq)  # original reads count

    # overlap reformat
    sfo = outdir + '/overlaps.sfo'
    overlaps = outdir + '/overlaps.savage'
    paf2sfo_cmd = 'python %s/minimap22sfo.py --in %s --out %s -m 0 -p 0' % (selfpath, paf, sfo)
    sfo2overlaps_cmd = 'python %s/sfo2overlaps.py --in %s --out %s --num_singles %d  --num_pairs 0' \
                       % (selfpath, sfo, overlaps, n_reads)
    subprocess.check_call(paf2sfo_cmd, shell=True)
    subprocess.check_call(sfo2overlaps_cmd, shell=True)
    # subprocess.check_call('rm -f %s' %sfo,shell=True)

    # construct digraph and simplify it
    run_vq = ' '.join([viralquasispecies,
                       "--singles", "%s " % single_fq,
                       "--overlaps=%s" % overlaps,
                       "--threads=%s" % threads,
                       "--graph_only=true",
                       "--edge_threshold=0",
                       "--ov_threshold=0",
                       "--min_qual=0",
                       "--cliques=true",
                       "--resolve_orientations=true",
                       "--error_correction=false",  # set flase to remove cycles
                       "--keep_singletons=0",
                       "--remove_branches=false",
                       "--remove_tips=%s" % rm_tips,
                       "--min_overlap_perc=0",
                       "--min_overlap_len=0",
                       "--FNO=3",
                       "--remove_trans=%d" % rm_trans,
                       "--optimize=false",
                       "--verbose=%s" % verbose,
                       "--base_path=%s" % (selfpath + '/../'),
                       "--min_read_len=%s" % min_read_len,
                       "--max_tip_len=%s" % max_tip_len,
                       "--separate_tips=true",
                       "--no_inclusion_overlaps=false",
                       "--ignore_inclusions=%s" % remove_inclusions,
                       "--original_readcount=%d" % n_reads,
                       "--output=%s/" % outdir
                       ]) + ' >/dev/null 2>&1 '
    subprocess.check_call(run_vq, shell=True)
    digraph_file = outdir + '/digraph.txt'
    return digraph_file


if __name__ == '__main__':
    fasta, paf, rm_trans, threads, remove_inclusions = sys.argv[1:]

    from assembly import rename_fa, rename_paf

    outdir = os.path.dirname(os.path.abspath(fasta))
    fasta2 = outdir + '/reads.renamed.fa'
    id_map_file = rename_fa(fasta, fasta2, outdir)

    paf2 = outdir + '/reads.renamed.paf'
    rename_paf(paf, paf2, id_map_file)
    ovlp2graph(fasta2, paf2, int(rm_trans), int(threads), remove_inclusions, 'true', 200, 2000)
