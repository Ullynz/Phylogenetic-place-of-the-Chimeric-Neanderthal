# Формирование химерного генома

from intervaltree import IntervalTree
import csv
import pandas as pd
import pysam

from config import PROJ_PATH, DAISEG_PATH, DAISEG_OUTGROUP
from support_functions import rnd, check_file

# Формируем химерный геном для современных популяций в участках неандертальской интрогрессии,
# используя гаплотип индивида с самым длинным архаичным участком, покрывающим очередной снип.
# Если снип не покрывается неандертальскими участками в данной популяции вообще, то в химерный геном 
# беру мажорный аллель
def form_modern_chimeric_genome(population, pop_A, pop_B):
    chimeric_genome_file = f"{PROJ_PATH}/{population}/{pop_A}.{pop_B}_chimeric_genome.tsv"

    if check_file(chimeric_genome_file):
        print(f"Chimeric genome for population {population} has already been built!")
        return
    
    with open(chimeric_genome_file, "w", newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["POP", "CHROM", "POS", "CHIMERIC_ALLELE", "REF_ALLELE"])
    
    for chrom in range(1, 23):
        tree = IntervalTree()
        daiseg_archaic_sections = f"{DAISEG_PATH}/{population}.{DAISEG_OUTGROUP}/{population}.{DAISEG_OUTGROUP}.grch37.chr{chrom}.em.tsv"
        df = pd.read_csv(daiseg_archaic_sections, sep='\t')

        for i in range(df.shape[0]):
            start = df.iloc[i]["Start"]
            end = df.iloc[i]["End"]
            sample = df.iloc[i]["Sample"]
            length = df.iloc[i]["Length"]
            tree.addi(start, end, (length, sample))

        vcf_name = f"{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"
        vcf_file = f"{PROJ_PATH}/{population}/{vcf_name}"

        vcf_data = pysam.VariantFile(vcf_file)
        with open(chimeric_genome_file, "a", newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            for rec in vcf_data:
                covering_sections = tree[rec.pos]
                if len(covering_sections):
                    longest_section = sorted(covering_sections, key=lambda x: x.data[0])[-1]
                    sample = longest_section.data[1]
                    allele = int(sample.split('_')[1]) - 1
                    sample = sample.split('_')[0]
                    if rec.samples[sample]['GT'][allele] == 1 and len(rec.alts[0]) == 1:
                        writer.writerow([population, rec.chrom, rec.pos, rec.alts[0], rec.ref])
                else:
                    AF = rec.info["AF"][0]
                    if AF > 0.5 and len(rec.alts[0]) == 1:
                        writer.writerow([population, rec.chrom, rec.pos, rec.alts[0], rec.ref])

# Формируем гаплоидный геном для неандертальцев 
def form_neand_chimeric_genome(population, pop_A, pop_B):
    chimeric_genome_file = f"{PROJ_PATH}/{population}/{pop_A}.{pop_B}_chimeric_genome.tsv"

    with open(chimeric_genome_file, "w", newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["POP", "CHROM", "POS", "CHIMERIC_ALLELE", "REF_ALLELE"])
    
    for chrom in range(1, 23):
        vcf_name = f"{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"
        vcf_file = f"{PROJ_PATH}/{population}/{vcf_name}"

        vcf_data = pysam.VariantFile(vcf_file)
        sample = vcf_data.header.samples[0]
        with open(chimeric_genome_file, "a", newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            for rec in vcf_data:
                allele = rnd()
                if rec.samples[sample]['GT'][allele] == 1 and len(rec.alts[0]) == 1:
                    writer.writerow([population, rec.chrom, rec.pos, rec.alts[0], rec.ref])

# Формирование "химерного генома" только для YRI (беру мажорные аллели), потом уберу этот блок после добавления генома Ust-ishim
def form_yri_chimeric_genome(population, pop_A, pop_B):
    chimeric_genome_file = f"{PROJ_PATH}/{population}/{pop_A}.{pop_B}_chimeric_genome.tsv"

    if check_file(chimeric_genome_file):
        print(f"Chimeric genome for population {population} has already been built!")
        return
    
    with open(chimeric_genome_file, "w", newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["POP", "CHROM", "POS", "CHIMERIC_ALLELE", "REF_ALLELE"])
    
    for chrom in range(1, 23):
        vcf_name = f"{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"
        vcf_file = f"{PROJ_PATH}/{population}/{vcf_name}"

        vcf_data = pysam.VariantFile(vcf_file)
        with open(chimeric_genome_file, "a", newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            for rec in vcf_data:
                AF = rec.info["AF"][0]
                if AF > 0.5 and len(rec.alts[0]) == 1:
                    writer.writerow([population, rec.chrom, rec.pos, rec.alts[0], rec.ref])