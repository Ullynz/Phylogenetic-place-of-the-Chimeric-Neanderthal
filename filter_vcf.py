# Фильтрация vcf по архаичным маскам

import subprocess

from config import PROJ_PATH, DAISEG_PATH, HUMAN_DATA, DAISEG_OUTGROUP
from support_functions import check_file, mkdir
from masks import merge_archaic_sections

ARCHAIC_SECTIONS_PATH = f"{PROJ_PATH}/merged_archaic_sections"

# Фильтр vcf современных популяций для конкретной хромосомы по сэмплам и маске
def modern_pop_filter(population, population_vcf, population_filtered_vcf, mask, samples, chrom):
    temp_vcf = f"{PROJ_PATH}/temp.vcf.gz"

    if not check_file(population_filtered_vcf):
        print(f"Filtering vcf for {population} chr{chrom}...")

        cmd = f"bcftools view -S {samples} -R {mask} --types snps {population_vcf} -Oz -o {temp_vcf}"
        subprocess.run(cmd, shell=True)

        cmd = f"bcftools +fill-tags {temp_vcf} -Oz -o {population_filtered_vcf}"
        subprocess.run(cmd, shell=True)

        cmd = f"bcftools index {population_filtered_vcf}"
        subprocess.run(cmd, shell=True)

        print("Done!")

# Фильтр vcf архаичных популяций для конкретной хромосомы по маске с отбрасыванием генотипа 0/0, который нам не интересен
def archaic_pop_filter(population, population_vcf, population_filtered_vcf, mask, chrom):
    if not check_file(population_filtered_vcf):
        print(f"Filtering vcf for {population} chr{chrom}...")

        cmd = f"""
        bcftools view -R {mask} --types snps {population_vcf} |
        bcftools view -e 'GT="0/0"' -Oz -o {population_filtered_vcf}
        """
        subprocess.run(cmd, shell=True)

        cmd = f"bcftools index {population_filtered_vcf}"
        subprocess.run(cmd, shell=True)

        print("Done!")

# Полная фильтрация
def filter_vcf(pop_A, pop_B):
    merge_archaic_sections(pop_A, pop_B)

    DAISEG_POP_A_PATH = f"{DAISEG_PATH}/{pop_A}.{DAISEG_OUTGROUP}"
    DAISEG_POP_B_PATH = f"{DAISEG_PATH}/{pop_B}.{DAISEG_OUTGROUP}"

    POP_A_PATH = f"{PROJ_PATH}/{pop_A}"
    POP_B_PATH = f"{PROJ_PATH}/{pop_B}"
    OUTGROUP_PATH = f"{PROJ_PATH}/{DAISEG_OUTGROUP}"
    VINDIJA_PATH = f"{PROJ_PATH}/Vindija"
    ALTAI_PATH = f"{PROJ_PATH}/Altai"

    mkdir(VINDIJA_PATH)
    mkdir(ALTAI_PATH)
    mkdir(OUTGROUP_PATH)

    # Извлекаем сэмплы из json для запуска DAIseg
    cmd = f"jq -r '.samples.ingroup[]' {DAISEG_POP_A_PATH}/jsons/{pop_A}.{DAISEG_OUTGROUP}.grch37.chr21.json > {POP_A_PATH}/samples.txt"
    subprocess.run(cmd, shell=True)
    pop_A_samples=f"{POP_A_PATH}/samples.txt"

    cmd = f"jq -r '.samples.ingroup[]' {DAISEG_POP_B_PATH}/jsons/{pop_B}.{DAISEG_OUTGROUP}.grch37.chr21.json > {POP_B_PATH}/samples.txt"
    subprocess.run(cmd, shell=True)
    pop_B_samples=f"{POP_B_PATH}/samples.txt"

    # Костыльно для YRI, пока строю с ним
    cmd = f"jq -r '.samples.outgroup[]' {DAISEG_POP_A_PATH}/jsons/{pop_A}.{DAISEG_OUTGROUP}.grch37.chr21.json > {OUTGROUP_PATH}/samples.txt"
    subprocess.run(cmd, shell=True)
    outgroup_samples=f"{OUTGROUP_PATH}/samples.txt"

    # Фильтруем
    for chrom in range(1, 23):
        Vindija_vcf = f"{HUMAN_DATA}/neand/33.19/chr{chrom}_mq25_mapab100.vcf.gz"
        Altai_vcf = f"{HUMAN_DATA}/neand/altai/chr{chrom}_mq25_mapab100.vcf.gz"
        pop_A_vcf = f"{DAISEG_POP_A_PATH}/1kG_filtered.chr{chrom}.grch37.vcf.gz"
        pop_B_vcf = f"{DAISEG_POP_B_PATH}/1kG_filtered.chr{chrom}.grch37.vcf.gz"
        outgroup_vcf = f"{HUMAN_DATA}/1000GP/1000GP.grch37/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

        Vindija_filtered_vcf = f"{VINDIJA_PATH}/{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"
        Altai_filtered_vcf = f"{ALTAI_PATH}/{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"
        pop_A_filtered_vcf = f"{POP_A_PATH}/{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"
        pop_B_filtered_vcf = f"{POP_B_PATH}/{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"
        outgroup_filtered_vcf = f"{OUTGROUP_PATH}/{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"

        mask = f"{ARCHAIC_SECTIONS_PATH}/{pop_A}.{pop_B}.chr{chrom}.bed"

        modern_pop_filter(pop_A, pop_A_vcf, pop_A_filtered_vcf, mask, pop_A_samples, chrom)
        modern_pop_filter(pop_B, pop_B_vcf, pop_B_filtered_vcf, mask, pop_B_samples, chrom)
        modern_pop_filter("YRI", outgroup_vcf, outgroup_filtered_vcf, mask, outgroup_samples, chrom)
        archaic_pop_filter("Vindija", Vindija_vcf, Vindija_filtered_vcf, mask, chrom)
        archaic_pop_filter("Altai", Altai_vcf, Altai_filtered_vcf, mask, chrom)
