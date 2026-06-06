# Модуль с вспомогательными функциями

from pathlib import Path
import random
import numpy as np
from collections import defaultdict

def rnd():
    return int(random.random() >= 0.5)

def divide(a, b):
    if b == 0:
        return 0
    return a / b

def list_to_string(lst):
    str = ".".join(obj for obj in lst)
    return str

# Работа с путями
def check_file(file_name):
    file_path = Path(file_name)

    if file_path.exists() and file_path.stat().st_size > 0:
        return True
    return False

def mkdir(dir_name):
    Path(dir_name).mkdir(parents=True, exist_ok=True)

# Удаление не биаллельных колонок из матрицы
def delete_bad_columns(matrix, populations):
    good_columns = []

    for pos in range(len(matrix[populations[0]])):
        variants = {}

        for population in populations:
            variants[matrix[population][pos]] = True

        if len(variants) == 2:
            good_columns.append(pos)
    
    result = {}
    for population in populations:
        result[population] = "".join(matrix[population][pos] for pos in good_columns)

    return result, good_columns

# Кодирование аллелей в 0/1
def code(allele, outgroup_allele):
    if allele == outgroup_allele:
        return 0
    return 1

# Статистики
def get_D_stat(matrix, population, start=0, end=None):
    if end is None:
        filtered_matrix, pos_idx = delete_bad_columns(matrix, [population, "Altai", "Vindija", "YRI"])
        end = len(filtered_matrix[population])

    else:
        filtered_matrix = matrix
        pos_idx = None

    end = min(end, len(filtered_matrix[population]))

    numer = 0.0
    denom = 0.0

    outgroup_seq = filtered_matrix["YRI"]
    for pos in range(start, end):

        outgroup_allele = outgroup_seq[pos]
        
        pop_coded = code(filtered_matrix[population][pos], outgroup_allele)
        vindija_coded = code(filtered_matrix["Vindija"][pos], outgroup_allele)
        altai_coded = code(filtered_matrix["Altai"][pos], outgroup_allele)

        value = -pop_coded * (altai_coded - vindija_coded)
        numer += value
        denom += abs(value)
    
    D_stat = divide(numer, denom)
    return numer, denom, D_stat, filtered_matrix, pos_idx

def get_Z_score_snp_blocks(matrix, population, numer, denom, D_stat):
    block_sz = 1250
    block_cnt = int((len(matrix[population]) - 1) / block_sz) + 1

    d_stats = []

    for block in range(block_cnt):
        start = block * block_sz
        end = (block + 1) * block_sz

        local_numer, local_denom, _, _, _ = get_D_stat(matrix, population, start, end)
        trunc_stat = divide(numer - local_numer, denom - local_denom)
        d_stats.append(trunc_stat)

    d_stats = np.array(d_stats)
    d_mean = np.mean(d_stats)

    SE = np.sqrt((block_cnt - 1) / block_cnt * np.sum((d_stats - d_mean)**2))
    Z_score = D_stat / SE

    return Z_score

def get_Z_score_physical_blocks(matrix, population, numer, denom, D_stat, pos, pos_idx):
    block_sz = 5e6
    blocks = defaultdict(list)

    for idx in pos_idx:
        blocks[(pos[idx][0], int(pos[idx][1] / block_sz))].append(idx)
    
    d_stats = []
    
    for block in blocks:
        start = blocks[block][0]
        end = blocks[block][-1] + 1
        
        local_numer, local_denom, _, _, _ = get_D_stat(matrix, population, start, end)
        trunc_stat = divide(numer - local_numer, denom - local_denom)
        d_stats.append(trunc_stat)
    
    d_stats = np.array(d_stats)
    d_mean = np.mean(d_stats)

    SE = np.sqrt((len(blocks) - 1) / len(blocks) * np.sum((d_stats - d_mean)**2))
    Z_score = D_stat / SE

    return Z_score