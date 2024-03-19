import sys
import os
from multiprocessing import Pool

from assembly import get_superead


def run_on_local(num_clusters, threads, outdir, type, min_cov, max_tip_len, n_correct, n_polish,
                 rm_trans, trim_ends,polish_tool,rm_tmp,correct_mode):
    params = []
    for j in range(num_clusters):
        i = str(j + 1)
        param = [i, outdir, type, min_cov, max_tip_len, n_correct, n_polish, rm_trans, trim_ends,polish_tool,rm_tmp,correct_mode]
        params.append(param)

    pool = Pool(threads)
    pool.map(get_superead, params, chunksize=1)  # ordered
    pool.close()
    pool.join()

    return


def run_on_hpc():
    ## TODO
    return
