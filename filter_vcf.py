# Фильтрация vcf по архаичным маскам

import subprocess
from multiprocessing import Pool

from config import PROJ_PATH, MODERN_POP, NEAND_POP, OUTGROUP, ARCHAIC_SECTIONS_PATH
from config import DAISEG_JSONS
from config import get_modern_pop_vcf, get_altai_vcf, get_vindija_vcf, get_outgroup_vcf
from support_functions import check_file, mkdir, list_to_string
from masks import intersect_archaic_sections

# Фильтр vcf современных популяций для конкретной хромосомы по сэмплам и маске
def modern_pop_filter(population, population_vcf, population_filtered_vcf, mask, samples, chrom):
    if check_file(population_filtered_vcf):
        return

    print(f"Filtering vcf for {population} chr{chrom}...")

    cmd = f"""
        bcftools view -R {mask} -S {samples} --types snps {population_vcf} |
        bcftools +fill-tags -Oz -o {population_filtered_vcf}
    """
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
def get_filtered_vcf(population_path, name, chrom):
    return f"{population_path}/{name}.filtered_chrom{chrom}.vcf.gz"

# Фильтрация всех популяций для заданной хромосомы
def filter_vcf_chrom(path, chrom):
    samples = {}
    vcf = {}
    filtered_vcf = {}

    for population in MODERN_POP + ["YRI"]:
        samples[population] = f"{path[population]}/samples.txt"
        vcf[population] = get_modern_pop_vcf(chrom)
    
    vcf["Vindija"] = get_vindija_vcf(chrom)
    vcf["Altai"] = get_altai_vcf(chrom)
    vcf[OUTGROUP] = get_outgroup_vcf(chrom)

    for population in MODERN_POP + NEAND_POP + [OUTGROUP] + ["YRI"]:
        filtered_vcf[population] = get_filtered_vcf(path[population], list_to_string(MODERN_POP), chrom)

    mask = f"{ARCHAIC_SECTIONS_PATH}/{list_to_string(MODERN_POP)}.chr{chrom}.bed"

    for population in MODERN_POP + ["YRI"]:
        modern_pop_filter(population, vcf[population], filtered_vcf[population], mask, samples[population], chrom)
   
    archaic_pop_filter("Vindija", vcf["Vindija"], filtered_vcf["Vindija"], mask, chrom)
    archaic_pop_filter("Altai", vcf["Altai"], filtered_vcf["Altai"], mask, chrom)

    if OUTGROUP != "YRI":
        archaic_pop_filter(OUTGROUP, vcf[OUTGROUP], filtered_vcf[OUTGROUP], mask, chrom)

# Полная фильтрация
def filter_vcf():
    intersect_archaic_sections()

    path = {}

    for population in MODERN_POP + NEAND_POP + [OUTGROUP] + ["YRI"]:
        path[population] = f"{PROJ_PATH}/{population}"
        mkdir(path[population])

    for idx, population in enumerate(MODERN_POP):
        pop_samples = f"{path[population]}/samples.txt"
        if not check_file(pop_samples):
            cmd = f"jq -r '.samples.ingroup[]' {DAISEG_JSONS[idx]} > {path[population]}/samples.txt"
            subprocess.run(cmd, shell=True)

    yri_samples=f"{path['YRI']}/samples.txt"
    if not check_file(yri_samples):
        cmd = f"jq -r '.samples.outgroup[]' {DAISEG_JSONS[0]} > {path['YRI']}/samples.txt"
        subprocess.run(cmd, shell=True)
    
    args = [(path, chrom) for chrom in range(1, 23)]

    with Pool(8) as pool:
        pool.starmap(filter_vcf_chrom, args)