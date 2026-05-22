# Работа с масками архаичных участков

import subprocess

from config import PROJ_PATH, DAISEG_PATH, DAISEG_OUTGROUP
from support_functions import mkdir, check_file

# Создание маски с архаичными участками + хорошим покрытием для популяции
def make_archaic_mask(population):
    DAISEG_POP_PATH = f"{DAISEG_PATH}/{population}.{DAISEG_OUTGROUP}"
    PROJ_POP_PATH = f"{PROJ_PATH}/{population}"

    mkdir(PROJ_POP_PATH)
    mkdir(f"{PROJ_POP_PATH}/archaic_masks")

    temp = f"{PROJ_POP_PATH}/archaic_masks/temp.bed"

    for chrom in range(1, 23):
        mask = f"{PROJ_POP_PATH}/archaic_masks/mask_chr{chrom}.bed"

        if not check_file(mask):
            print(f"Forming mask for {population} chrom {chrom}...")

            coverage = f"{DAISEG_POP_PATH}/coverage_1kG.chr{chrom}.grch37.bed"
            good_coverage_mask = f"{PROJ_POP_PATH}/archaic_masks/coverage_mask_chr{chrom}.bed"
            cmd = f'awk \'$7 > 0.8 {{print $1"\t"$2"\t"$3}}\' {coverage} > {good_coverage_mask}'
            subprocess.run(cmd, shell=True)

            archaic_sections = f"{DAISEG_POP_PATH}/{population}.{DAISEG_OUTGROUP}.grch37.chr{chrom}.em.tsv"
            archaic_sections_mask = f"{PROJ_POP_PATH}/archaic_masks/archaic_mask_chr{chrom}.bed"
            cmd = f'awk \'{{print $2"\t"$3"\t"$4}}\' {archaic_sections} > {archaic_sections_mask}'
            subprocess.run(cmd, shell=True)

            cmd = f"bedtools intersect -a {good_coverage_mask} -b {archaic_sections_mask} > {temp}"
            subprocess.run(cmd, shell=True)

            cmd = f"sort -k1,1 -k2,2n {temp} | bedtools merge > {mask}"
            subprocess.run(cmd, shell=True)

            print("Done!")

# Объединение архаичных участков для выбранных популяций
def merge_archaic_sections(pop_A, pop_B):
    make_archaic_mask(pop_A)
    make_archaic_mask(pop_B)

    PROJ_POP_A_PATH=f"{PROJ_PATH}/{pop_A}"
    PROJ_POP_B_PATH=f"{PROJ_PATH}/{pop_B}"
    ARCHAIC_SECTIONS_PATH=f"{PROJ_PATH}/merged_archaic_sections"

    mkdir(ARCHAIC_SECTIONS_PATH)

    temp=f"{ARCHAIC_SECTIONS_PATH}/temp.bed"

    for chrom in range(1, 23):
        archaic_sections_A=f"{PROJ_POP_A_PATH}/archaic_masks/mask_chr{chrom}.bed"
        archaic_sections_B=f"{PROJ_POP_B_PATH}/archaic_masks/mask_chr{chrom}.bed"
        merged_sections=f"{ARCHAIC_SECTIONS_PATH}/{pop_A}.{pop_B}.chr{chrom}.bed"

        if not check_file(merged_sections):
            print(f"Merging archaic sections of {pop_A} and {pop_B} chr{chrom}...")

            cmd = f"cat {archaic_sections_A} > {temp}"
            subprocess.run(cmd, shell=True)

            cmd = f"cat {archaic_sections_B} >> {temp}"
            subprocess.run(cmd, shell=True)
            
            cmd = f"sort -k1,1 -k2,2n {temp} | bedtools merge > {merged_sections}"
            subprocess.run(cmd, shell=True)

            print("Done!")    