# Построение MSA матрицы

import sys
import pandas as pd
from collections import defaultdict

from config import PROJ_PATH
from support_functions import check_file

MSA_MATRIX_FASTA = f"{PROJ_PATH}/msa_matrix.fasta"

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

# Строим MSA матрицу. Места пропусков (т.е. места отсутствия снипов в какой-то из популяций на данной позиции)
# заполяняю референсным алеллем.
def form_msa_matrix(populations):
    pop_A = populations[0]
    pop_B = populations[1]

    matrix = defaultdict(dict)
    ordered_matrix = defaultdict(dict)
    ordered_matrix_str = {}

    # Заполняю всю матрицу референсными аллелями
    for population in populations:
        chimeric_genome = f"{PROJ_PATH}/{population}/{pop_A}.{pop_B}_chimeric_genome.tsv"

        if not check_file(chimeric_genome):
            print(f"Chimeric genome for population {population} was not found")
            sys.exit(1)

        df = pd.read_csv(chimeric_genome, sep='\t')

        for i in range(df.shape[0]):
            ref = df.iloc[i]["REF_ALLELE"]
            chrom = df.iloc[i]["CHROM"]
            pos = df.iloc[i]["POS"]

            for population_ in populations:
                matrix[population_][(chrom, pos)] = ref
    
    # Заменяю в матрице референсные аллели на снипы, где они есть
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

    # Перевожу в fasta формат
    for population in populations:
        ordered_matrix_str[population] = "".join(ordered_matrix[population].values())
    
    out_fasta = MSA_MATRIX_FASTA
    with open(out_fasta, "w") as f:
        for population in populations:
            f.write(f">{population} chimeric genome\n")
            f.write(f"{ordered_matrix_str[population]}\n")