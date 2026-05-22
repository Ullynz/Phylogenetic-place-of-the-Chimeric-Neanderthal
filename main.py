import sys

from config import PROJ_PATH, NEAND_POP, TREE_OUTGROUP
from filter_vcf import filter_vcf
from tree import run_iqtree
from msa_matrix import form_msa_matrix
from chimeric_genome import form_modern_chimeric_genome, form_neand_chimeric_genome, form_yri_chimeric_genome

MSA_MATRIX_FASTA = f"{PROJ_PATH}/msa_matrix.fasta"

# Главная функция, ожидает на вход 2 современные популяции, которые будем рассматривать
def main():
    if len(sys.argv) != 3:
        print("ERROR: 2 populations should be given")
        sys.exit(1)

    modern_populations = sys.argv[1:]
    modern_populations.sort(reverse=True)

    pop_A = modern_populations[0]
    pop_B = modern_populations[1]

    filter_vcf(pop_A, pop_B)

    for population in modern_populations:
        print(f"Processing of population {population}...")
        form_modern_chimeric_genome(population, pop_A, pop_B)
        print("Chimeric genome was formed successfully!")
    
    for population in NEAND_POP:        
        print(f"Processing of population {population}...")
        form_neand_chimeric_genome(population, pop_A, pop_B)
        print("Chimeric genome was formed successfully!")
    
    # Временный блок с YRI, потом удалю его после добавления Ust-ishim
    print(f"Processing of population YRI...")
    form_yri_chimeric_genome("YRI", pop_A, pop_B)
    print("Chimeric genome was formed successfully!")

    print("Forming MSA matrix...")
    all_populations = modern_populations + NEAND_POP + [TREE_OUTGROUP]
    form_msa_matrix(all_populations)
    print("Done!")

    print("Building a tree...")
    run_iqtree(MSA_MATRIX_FASTA)

    print("All jobs were finished!")

if __name__ == "__main__":
    main()
    