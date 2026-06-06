from config import MODERN_POP, NEAND_POP, OUTGROUP, DAISEG_PATHS
from filter_vcf import filter_vcf
from tree import run_iqtree
from msa_matrix import form_msa_matrix
from chimeric_genome import form_modern_chimeric_genome, form_neand_chimeric_genome, form_yri_genome
from support_functions import get_D_stat, get_Z_score_snp_blocks, get_Z_score_physical_blocks

def main():
    filter_vcf()

    print("Forming chimeric genome:\n")

    for idx, population in enumerate(MODERN_POP):
        print(f"Processing of population {population}...")
        form_modern_chimeric_genome(population, DAISEG_PATHS[idx])
   
    for population in NEAND_POP + [OUTGROUP]:
        if population == "YRI":
            continue

        print(f"Processing of population {population}...")
        form_neand_chimeric_genome(population)
    
    print(f"Processing of population YRI...")
    form_yri_genome()

    print("Forming MSA matrix...")
    all_populations = MODERN_POP + NEAND_POP + [OUTGROUP]
    matrix, pos = form_msa_matrix(all_populations)
    print("Done!\n")

    print("Calculating statistics...")
    for population in MODERN_POP:
        numer, denom, D_stat, filtered_matrix, pos_idx = get_D_stat(matrix, population)
        print(f"D(YRI, {population}; Altai, Vindija) = {D_stat}")

        Z_score_snp = get_Z_score_snp_blocks(filtered_matrix, population, numer, denom, D_stat)
        print(f"Z_snp(YRI, {population}; Altai, Vindija) = {Z_score_snp}")

        Z_score_physical = get_Z_score_physical_blocks(filtered_matrix, population, numer, denom, D_stat, pos, pos_idx)
        print(f"Z_physical(YRI, {population}; Altai, Vindija) = {Z_score_physical}")
    print()

    print("Building trees...")
    run_iqtree()

    print("All jobs were finished!")

if __name__ == "__main__":
    main()
    