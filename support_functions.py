# Модуль с вспомогательными функциями

from pathlib import Path
import random
import numpy as np

def rnd():
    return int(random.random() >= 0.5)

def check_file(file_name):
    file_path = Path(file_name)

    if file_path.exists() and file_path.stat().st_size > 0:
        return True
    return False

def mkdir(dir_name):
    Path(dir_name).mkdir(parents=True, exist_ok=True)

def code(allele, outgroup_allele):
    if allele == outgroup_allele:
        return 0
    return 1

def get_D_stat(matrix, population, outgroup, start=0, end=None):
    if end is None:
        end = len(matrix[population])
    
    numer = 0.0
    denom = 0.0
    outgroup_seq = matrix[outgroup]
    for pos in range(start, end):
        outgroup_allele = outgroup_seq[pos]
        
        pop_coded = code(matrix[population][pos], outgroup_allele)
        vindija_coded = code(matrix["Vindija"][pos], outgroup_allele)
        altai_coded = code(matrix["Altai"][pos], outgroup_allele)

        value = -pop_coded * (altai_coded - vindija_coded)
        numer += value
        denom += abs(value)
    
    D_stat = numer / denom
    return numer, denom, D_stat

def get_Z_score(matrix, population, outgroup, numer, denom, D_stat):
    block_sz = 5000
    block_cnt = int((len(matrix[population]) - 1) / block_sz)

    d_stats = []

    for block in range(block_cnt):
        start = block * block_sz
        end = min((block + 1) * block_sz, len(matrix[population]))

        local_numer, local_denom, _ = get_D_stat(matrix, population, outgroup, start, end)
        trunc_stat = (numer - local_numer) / (denom - local_denom)
        d_stats.append(trunc_stat)

    d_stats = np.array(d_stats)
    d_mean = np.mean(d_stats)

    SE = np.sqrt((block_cnt - 1) / block_cnt * np.sum((d_stats - d_mean)**2))
    Z_score = D_stat / SE

    return Z_score

