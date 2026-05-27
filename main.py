from config import MODERN_POP, NEAND_POP, OUTGROUP
from config import get_daiseg_pop_A_output, get_daiseg_pop_B_output
from filter_vcf import filter_vcf
from tree import run_iqtree
from msa_matrix import form_msa_matrix
from chimeric_genome import form_modern_chimeric_genome, form_neand_chimeric_genome

def main():
    pop_A = MODERN_POP[0]
    pop_B = MODERN_POP[1]

    filter_vcf(pop_A, pop_B)

    print("Forming chimeric genome:\n")

    print(f"Processing of population {pop_A}...")
    form_modern_chimeric_genome(population=pop_A, pop_A=pop_A, pop_B=pop_B, get_daiseg_output=get_daiseg_pop_A_output)

    print(f"Processing of population {pop_B}...")
    form_modern_chimeric_genome(population=pop_B, pop_A=pop_A, pop_B=pop_B, get_daiseg_output=get_daiseg_pop_B_output)
    
    for population in NEAND_POP + [OUTGROUP]:        
        print(f"Processing of population {population}...")
        form_neand_chimeric_genome(population=population, pop_A=pop_A, pop_B=pop_B)

    print("Forming MSA matrix...")
    all_populations = MODERN_POP + NEAND_POP + [OUTGROUP]
    form_msa_matrix(all_populations)
    print("Done!\n")

    print("Building trees...")
    run_iqtree()

    print("All jobs were finished!")

if __name__ == "__main__":
    main()
    