# Работа с масками архаичных участков

import subprocess

from config import PROJ_PATH, DAISEG_POP_A_PATH, DAISEG_POP_B_PATH, ARCHAIC_SECTIONS_PATH
from config import get_daiseg_coverage, get_daiseg_pop_A_output, get_daiseg_pop_B_output, get_altai_coverage, get_vindija_coverage, get_outgroup_coverage
from support_functions import mkdir, check_file

# Создание маски с архаичными участками + хорошим покрытием для популяции
def make_archaic_mask(population, daiseg_path, get_daiseg_output):
    PROJ_POP_PATH = f"{PROJ_PATH}/{population}"

    mkdir(PROJ_POP_PATH)
    mkdir(f"{PROJ_POP_PATH}/archaic_masks")

    temp = f"{PROJ_POP_PATH}/archaic_masks/temp.bed"

    for chrom in range(1, 23):
        mask = f"{PROJ_POP_PATH}/archaic_masks/mask_chr{chrom}.bed"

        if check_file(mask):
            continue

        print(f"Forming mask for {population} chrom {chrom}...")

        coverage = get_daiseg_coverage(daiseg_path, chrom)
        good_coverage_mask = f"{PROJ_POP_PATH}/archaic_masks/coverage_mask_chr{chrom}.bed"
        cmd = f'awk \'$7 > 0.8 {{print $1"\t"$2"\t"$3}}\' {coverage} > {good_coverage_mask}'
        subprocess.run(cmd, shell=True)

        archaic_sections = get_daiseg_output(chrom)
        archaic_sections_mask = f"{PROJ_POP_PATH}/archaic_masks/archaic_mask_chr{chrom}.bed"
        cmd = f'awk \'{{print $2"\t"$3"\t"$4}}\' {archaic_sections} > {archaic_sections_mask}'
        subprocess.run(cmd, shell=True)

        cmd = f"bedtools intersect -a {good_coverage_mask} -b {archaic_sections_mask} > {temp}"
        subprocess.run(cmd, shell=True)

        cmd = f"sort -k1,1 -k2,2n {temp} | bedtools merge > {mask}"
        subprocess.run(cmd, shell=True)

        print("Done!")

# Пересение масок покрытия неандретальских геномов
def intersect_neand_masks():
    temp=f"{ARCHAIC_SECTIONS_PATH}/temp.bed"

    for chrom in range(1, 23):
        altai_mask = get_altai_coverage(chrom)
        vindija_mask = get_vindija_coverage(chrom)
        intersected_masks = f"{ARCHAIC_SECTIONS_PATH}/neand_mask_chr{chrom}.bed"

        if check_file(intersected_masks):
            continue
        
        print(f"Intersecting neand coverage masks for chr{chrom}...")

        cmd = f"bedtools intersect -a {altai_mask} -b {vindija_mask} > {temp}"
        subprocess.run(cmd, shell=True)

        cmd = f"sort -k1,1 -k2,2n {temp} | bedtools merge > {intersected_masks}"
        subprocess.run(cmd, shell=True)

        print("Done!")

# Пересечение масок всех 5 рассматриваемых популяций
def intersect_archaic_sections(pop_A, pop_B):
    make_archaic_mask(pop_A, DAISEG_POP_A_PATH, get_daiseg_pop_A_output)
    make_archaic_mask(pop_B, DAISEG_POP_B_PATH, get_daiseg_pop_B_output)
    intersect_neand_masks()

    PROJ_POP_A_PATH = f"{PROJ_PATH}/{pop_A}"
    PROJ_POP_B_PATH = f"{PROJ_PATH}/{pop_B}"

    mkdir(ARCHAIC_SECTIONS_PATH)

    temp=f"{ARCHAIC_SECTIONS_PATH}/temp.bed"
    temp2=f"{ARCHAIC_SECTIONS_PATH}/temp2.bed"

    for chrom in range(1, 23):
        mask_A = f"{PROJ_POP_A_PATH}/archaic_masks/mask_chr{chrom}.bed"
        mask_B = f"{PROJ_POP_B_PATH}/archaic_masks/mask_chr{chrom}.bed"
        neand_mask = f"{ARCHAIC_SECTIONS_PATH}/neand_mask_chr{chrom}.bed"
        outgroup_mask = get_outgroup_coverage(chrom)

        intersected_sections = f"{ARCHAIC_SECTIONS_PATH}/{pop_A}.{pop_B}.chr{chrom}.bed"

        if check_file(intersected_sections):
            continue
        
        print(f"Intersecting archaic sections of {pop_A} and {pop_B} chr{chrom}...")

        cmd = f"bedtools intersect -a {mask_A} -b {mask_B} > {temp}"
        subprocess.run(cmd, shell=True)

        cmd = f"bedtools intersect -a {temp} -b {neand_mask} > {temp2}"
        subprocess.run(cmd, shell=True)

        cmd = f"bedtools intersect -a {temp2} -b {outgroup_mask} > {temp}"
        subprocess.run(cmd, shell=True)

        cmd = f"sort -k1,1 -k2,2n {temp} | bedtools merge > {intersected_sections}"
        subprocess.run(cmd, shell=True)

        print("Done!")