# Фильтрация vcf по архаичным маскам

import subprocess

from config import PROJ_PATH, OUTGROUP, ARCHAIC_SECTIONS_PATH
from config import DAISEG_POP_A_JSON, DAISEG_POP_B_JSON, DAISEG_POP_A_PATH, DAISEG_POP_B_PATH
from config import get_modern_pop_vcf, get_altai_vcf, get_vindija_vcf, get_outgroup_vcf, get_yri_vcf
from support_functions import check_file, mkdir
from masks import intersect_archaic_sections

# Фильтр vcf современных популяций для конкретной хромосомы по сэмплам и маске
def modern_pop_filter(population, population_vcf, population_filtered_vcf, mask, samples, chrom):
    temp_vcf = f"{PROJ_PATH}/temp.vcf.gz"

    if check_file(population_filtered_vcf):
        return

    print(f"Filtering vcf for {population} chr{chrom}...")

    cmd = f"bcftools view -R {mask} -S {samples} --types snps {population_vcf} -Oz -o {temp_vcf}"
    subprocess.run(cmd, shell=True)

    cmd = f"bcftools +fill-tags {temp_vcf} -Oz -o {population_filtered_vcf}"
    subprocess.run(cmd, shell=True)

    cmd = f"bcftools index {population_filtered_vcf}"
    subprocess.run(cmd, shell=True)

    print("Done!")

# Фильтр vcf архаичных популяций для конкретной хромосомы по маске сразу с отбрасыванием генотипа 0/0
def archaic_pop_filter(population, population_vcf, population_filtered_vcf, mask, chrom):
    if check_file(population_filtered_vcf):
        return
    
    print(f"Filtering vcf for {population} chr{chrom}...")

    cmd = f"bcftools view -R {mask} --types snps -e 'GT=\"0/0\"' {population_vcf} -Oz -o {population_filtered_vcf}"
    subprocess.run(cmd, shell=True)

    cmd = f"bcftools index {population_filtered_vcf}"
    subprocess.run(cmd, shell=True)

    print("Done!")

# Путь к отфильтрованным по архаичным сегментам vcf
def get_filtered_vcf(population_path, pop_A, pop_B, chrom):
    return f"{population_path}/{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"

# Полная фильтрация
def filter_vcf(pop_A, pop_B):
    intersect_archaic_sections(pop_A, pop_B)

    POP_A_PATH = f"{PROJ_PATH}/{pop_A}"
    POP_B_PATH = f"{PROJ_PATH}/{pop_B}"
    OUTGROUP_PATH = f"{PROJ_PATH}/{OUTGROUP}"
    VINDIJA_PATH = f"{PROJ_PATH}/Vindija"
    ALTAI_PATH = f"{PROJ_PATH}/Altai"
    YRI_PATH = f"{PROJ_PATH}/YRI"

    mkdir(VINDIJA_PATH)
    mkdir(ALTAI_PATH)
    mkdir(OUTGROUP_PATH)
    mkdir(YRI_PATH)

    cmd = f"jq -r '.samples.ingroup[]' {DAISEG_POP_A_JSON} > {POP_A_PATH}/samples.txt"
    subprocess.run(cmd, shell=True)
    pop_A_samples=f"{POP_A_PATH}/samples.txt"

    cmd = f"jq -r '.samples.ingroup[]' {DAISEG_POP_B_JSON} > {POP_B_PATH}/samples.txt"
    subprocess.run(cmd, shell=True)
    pop_B_samples=f"{POP_B_PATH}/samples.txt"

    cmd = f"jq -r '.samples.outgroup[]' {DAISEG_POP_B_JSON} > {YRI_PATH}/samples.txt"
    subprocess.run(cmd, shell=True)
    yri_samples=f"{YRI_PATH}/samples.txt"

    for chrom in range(1, 23):
        Vindija_vcf = get_vindija_vcf(chrom)
        Altai_vcf = get_altai_vcf(chrom)
        pop_A_vcf = get_modern_pop_vcf(DAISEG_POP_A_PATH, chrom)
        pop_B_vcf = get_modern_pop_vcf(DAISEG_POP_B_PATH, chrom)
        outgroup_vcf = get_outgroup_vcf(chrom)
        yri_vcf = get_yri_vcf(chrom)

        Vindija_filtered_vcf = get_filtered_vcf(VINDIJA_PATH, pop_A, pop_B, chrom)
        Altai_filtered_vcf = get_filtered_vcf(ALTAI_PATH, pop_A, pop_B, chrom)
        pop_A_filtered_vcf = get_filtered_vcf(POP_A_PATH, pop_A, pop_B, chrom)
        pop_B_filtered_vcf = get_filtered_vcf(POP_B_PATH, pop_A, pop_B, chrom)
        outgroup_filtered_vcf = get_filtered_vcf(OUTGROUP_PATH, pop_A, pop_B, chrom)
        yri_filtered_vcf = get_filtered_vcf(YRI_PATH, pop_A, pop_B, chrom)

        mask = f"{ARCHAIC_SECTIONS_PATH}/{pop_A}.{pop_B}.chr{chrom}.bed"

        modern_pop_filter(pop_A, pop_A_vcf, pop_A_filtered_vcf, mask, pop_A_samples, chrom)
        modern_pop_filter(pop_B, pop_B_vcf, pop_B_filtered_vcf, mask, pop_B_samples, chrom)
        archaic_pop_filter("Vindija", Vindija_vcf, Vindija_filtered_vcf, mask, chrom)
        archaic_pop_filter("Altai", Altai_vcf, Altai_filtered_vcf, mask, chrom)
        archaic_pop_filter(OUTGROUP, outgroup_vcf, outgroup_filtered_vcf, mask, chrom)
        modern_pop_filter("YRI", yri_vcf, yri_filtered_vcf, mask, yri_samples, chrom)