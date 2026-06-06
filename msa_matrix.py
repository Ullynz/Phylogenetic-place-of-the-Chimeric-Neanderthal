# Построение MSA матрицы

import pandas as pd
from collections import defaultdict

from config import PROJ_PATH, MODERN_POP, MSA_MATRIX_PATH, TREES_NUM, WINDOW_SIZE
from support_functions import mkdir, delete_bad_columns, list_to_string

def form_msa_matrix(populations):
    matrix = defaultdict(dict)
    ordered_matrix = defaultdict(dict)
    ordered_matrix_str = {}
    filtered_ordered_matrix_str = {}

    for population in populations + ["YRI"]:
        chimeric_genome = f"{PROJ_PATH}/{population}/{list_to_string(MODERN_POP)}_chimeric_genome.tsv"
        df = pd.read_csv(chimeric_genome, sep='\t')

        for i in range(df.shape[0]):
            ref = df.iloc[i]["REF_ALLELE"]
            chrom = df.iloc[i]["CHROM"]
            pos = df.iloc[i]["POS"]

            for population_ in populations + ["YRI"]:
                matrix[population_][(chrom, pos)] = ref
    
    for population in populations + ["YRI"]:
        chimeric_genome = f"{PROJ_PATH}/{population}/{list_to_string(MODERN_POP)}_chimeric_genome.tsv"
        df = pd.read_csv(chimeric_genome, sep='\t')

        for i in range(df.shape[0]):
            chimeric_allele = df.iloc[i]["CHIMERIC_ALLELE"]
            chrom = df.iloc[i]["CHROM"]
            pos = df.iloc[i]["POS"]

            matrix[population][(chrom, pos)] = chimeric_allele

        ordered_matrix[population] = dict(sorted(matrix[population].items()))

    for population in populations + ["YRI"]:
        ordered_matrix_str[population] = "".join(ordered_matrix[population].values())
    
    filtered_ordered_matrix_str, _ = delete_bad_columns(ordered_matrix_str, populations)

    mkdir(MSA_MATRIX_PATH)

    matrix_len = len(filtered_ordered_matrix_str[populations[0]])
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

                f.write(f"{filtered_ordered_matrix_str[population][start:end]}\n")
    
    return ordered_matrix_str, list(ordered_matrix[populations[0]].keys())