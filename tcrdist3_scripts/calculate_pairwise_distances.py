#!/usr/bin/python3

import pandas as pd
from tcrdist.repertoire import TCRrep
from tcrdist.adpt_funcs import import_adaptive_file, adaptive_to_imgt
import math
import os
import time
from scipy.sparse import triu

def get_safe_chunk(clones, cpus = 40, target = 7*(10**7)):
    max_oper_sim = target * cpus
    tot_oper = clones**2

    if tot_oper < max_oper_sim:
        target = tot_oper / cpus

    ideal_chunk_size = math.ceil(target / clones)
    return ideal_chunk_size



target = 7*(10**7)
cpus = 40



os.chdir('janssen_vaccine/')
files = os.listdir(os.getcwd())
files.sort(key=lambda f: os.stat(f).st_size, reverse=False)
files = [f for f in files if f.endswith('.tsv')]
# files = files[0:2] ############# TEST : smallest files
# files = [files[len(files)-1]] ############# TEST : largest file (max RAM 215GB)
# print(files)
tot_start = time.time()

for f in files:
    log = open("../output/run_tcrdist3.log", "a")

    log.write("\n#####################################\n" + f + "\n#####################################")
    start = time.time()

    df = import_adaptive_file(adaptive_filename = f, # 'janssen_vaccine/Vaccine_Subj01.tsv' (Largest sample !!!!)
                          organism = "human",
                          chain = "beta",
                          count = 'productive_frequency')

    tr = TCRrep(cell_df = df,
                organism = 'human',
                chains = ['beta'],
                db_file = 'alphabeta_gammadelta_db.tsv',
                compute_distances = False,  # For larger datasets, make sure compute_distances is set to False,
                store_all_cdr = False)

    tr.cpus = cpus

    log.write("\n### Total clones: " + str(tr.clone_df.shape[0]))

    safe_chunk_size = get_safe_chunk(
                tr.clone_df.shape[0],
                cpus = cpus,
                target = target)

    log.write("\n### Chunk size: " + str(safe_chunk_size))


    tr.compute_sparse_rect_distances(
            df = tr.clone_df,
            radius=12,
            chunk_size = safe_chunk_size)

    csc = triu(tr.rw_beta) # Upper triangular
    coo = csc.tocoo(copy=False)
    df = pd.DataFrame({'from_idx': coo.row, 'to_idx': coo.col, 'distance': coo.data}
                     )[['from_idx', 'to_idx', 'distance']].sort_values(['from_idx', 'to_idx']
                     ).reset_index(drop=True)

    # df = df[df.distance != -1]
    df = df[df['from_idx'] != df['to_idx']] # Remove pw dist between the same clone

    clone_info = tr.clone_df[['cdr3_b_aa', 'v_b_gene', 'templates']]
    clone_info['idx'] = clone_info.index


    # df['from_cdr3'] = clone_info.loc[list(df['from_idx'])].cdr3_b_aa.tolist()
    # df['to_cdr3'] = clone_info.loc[list(df['to_idx'])].cdr3_b_aa.tolist()
    # df['from_v'] = clone_info.loc[list(df['from_idx'])].v_b_gene.tolist()
    # df['to_v'] = clone_info.loc[list(df['to_idx'])].v_b_gene.tolist()

    df.to_csv('../output/' + os.path.splitext(f)[0] + '_dist.csv', index = False)
    clone_info.to_csv('../output/' + os.path.splitext(f)[0] + '_cloneinfo.csv', index = False)

    end = time.time()
    log.write("\n### Runtime: " + str(end - start))
    log.close()

tot_end = time.time()
log = open("../output/run_tcrdist3.log", "a")
log.write("\n######################################\n### TOTAL Runtime: " + str(tot_end - tot_start))
log.close()
