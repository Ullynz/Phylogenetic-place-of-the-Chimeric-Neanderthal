# Формирование химерного генома

from intervaltree import IntervalTree
import csv
import pandas as pd
import pysam

from config import PROJ_PATH
from support_functions import rnd, check_file

# Формируем химерный геном для современных популяций в участках неандертальской интрогрессии,
# используя гаплотип индивида с самым длинным архаичным участком, покрывающим очередной снип.
def form_modern_chimeric_genome(population, pop_A, pop_B, get_daiseg_output):
    chimeric_genome_file = f"{PROJ_PATH}/{population}/{pop_A}.{pop_B}_chimeric_genome.tsv"

    if check_file(chimeric_genome_file):
        print(f"Chimeric genome for population {population} has already been built!\n")
        return

    neand_snp = 0
    non_neand_snp = 0
    
    with open(chimeric_genome_file, "w", newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["POP", "CHROM", "POS", "CHIMERIC_ALLELE", "REF_ALLELE"])
    
    for chrom in range(1, 23):
        tree = IntervalTree()
        daiseg_archaic_sections = get_daiseg_output(chrom)
        df = pd.read_csv(daiseg_archaic_sections, sep='\t')

        for i in range(df.shape[0]):
            start = df.iloc[i]["Start"]
            end = df.iloc[i]["End"]
            sample = df.iloc[i]["Sample"]
            length = df.iloc[i]["Length"]
            tree.addi(start, end + 1, (length, sample))

        vcf_name = f"{pop_A}.{pop_B}.filtered_chrom{chrom}.vcf.gz"
        vcf_file = f"{PROJ_PATH}/{population}/{vcf_name}"

        vcf_data = pysam.VariantFile(vcf_file)
        with open(chimeric_genome_file, "a", newline='') as f:
            writer = csv.writer(f, delimiter='\t')

            for rec in vcf_data:
                covering_sections = tree[rec.pos]

                if len(covering_sections) > 0:
                    longest_section = sorted(covering_sections, key=lambda x: x.data[0])[-1]
                    sample = longest_section.data[1]
                    allele = int(sample.split('_')[1]) - 1
                    sample = sample.split('_')[0]

                    if rec.samples[sample]['GT'][allele] == 1 and len(rec.alts[0]) == 1:
                        neand_snp += 1
                        writer.writerow([population, rec.chrom, rec.pos, rec.alts[0], rec.ref])
                else:
                    AF = rec.info["AF"][0]
                    if AF > 0.5 and len(rec.alts[0]) == 1:
                        non_neand_snp += 1
                        writer.writerow([population, rec.chrom, rec.pos, rec.alts[0], rec.ref])

    total_snp = neand_snp + non_neand_snp
    neand_snp_percentage = neand_snp * 100 / (total_snp)

    print(f"~{neand_snp_percentage}% of {population} snp are located in Neanderthal sections, {non_neand_snp} major alleles were used")
    print(f"Chimeric genome was formed successfully! {total_snp} snp were added\n")

# Формирование гаплоидного генома для неандертальцев и outgroup
def form_neand_chimeric_genome(population, pop_A, pop_B):
    chimeric_genome_file = f"{PROJ_PATH}/{population}/{pop_A}.{pop_B}_chimeric_genome.tsv"

    total_snp = 0

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
                    total_snp += 1
                    writer.writerow([population, rec.chrom, rec.pos, rec.alts[0], rec.ref])
    
    print(f"Chimeric genome was formed successfully! {total_snp} snp were added\n")