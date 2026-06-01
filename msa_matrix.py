# Построение MSA матрицы

import pandas as pd
from collections import defaultdict

from config import PROJ_PATH, MSA_MATRIX_PATH, TREES_NUM, WINDOW_SIZE
from support_functions import mkdir

# Удаление не биаллельных колонок из матрицы
def delete_bad_columns(matrix, populations):
    bad_columns = []

    for pos in matrix[populations[0]]:
        variants = {}

        for population in populations:
            variants[matrix[population][pos]] = True

        if len(variants) != 2:
            bad_columns.append(pos)
    
    for population in populations:
        for col in bad_columns:
            del matrix[population][col]

    return matrix

# Построение MSA матрицы
def form_msa_matrix(populations):
    pop_A = populations[0]
    pop_B = populations[1]

    matrix = defaultdict(dict)
    ordered_matrix = defaultdict(dict)
    ordered_matrix_str = {}

    for population in populations:
        chimeric_genome = f"{PROJ_PATH}/{population}/{pop_A}.{pop_B}_chimeric_genome.tsv"
        df = pd.read_csv(chimeric_genome, sep='\t')

        for i in range(df.shape[0]):
            ref = df.iloc[i]["REF_ALLELE"]
            chrom = df.iloc[i]["CHROM"]
            pos = df.iloc[i]["POS"]

            for population_ in populations:
                matrix[population_][(chrom, pos)] = ref
    
    for population in populations:
        chimeric_genome = f"{PROJ_PATH}/{population}/{pop_A}.{pop_B}_chimeric_genome.tsv"
        df = pd.read_csv(chimeric_genome, sep='\t')

        for i in range(df.shape[0]):
            chimeric_allele = df.iloc[i]["CHIMERIC_ALLELE"]
            chrom = df.iloc[i]["CHROM"]
            pos = df.iloc[i]["POS"]

            matrix[population][(chrom, pos)] = chimeric_allele

        ordered_matrix[population] = dict(sorted(matrix[population].items()))
    
    ordered_matrix = delete_bad_columns(ordered_matrix, populations)

    for population in populations:
        ordered_matrix_str[population] = "".join(ordered_matrix[population].values())

    mkdir(MSA_MATRIX_PATH)

    matrix_len = len(ordered_matrix_str[populations[0]])
    block_size = int(matrix_len / TREES_NUM) + 1
    
    for cnt in range(TREES_NUM):
        out_fasta = f"{MSA_MATRIX_PATH}/matrix_{cnt}.fasta"
        with open(out_fasta, "w") as f:
            for population in populations:
                f.write(f">{population} chimeric genome\n")

                start = cnt * block_size
                end = min((cnt + 1) * block_size, matrix_len)
                
                if WINDOW_SIZE > 0:
                    end = min(end, start + WINDOW_SIZE)

                f.write(f"{ordered_matrix_str[population][start:end]}\n")
    
    return ordered_matrix_str