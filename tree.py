# Строительство филогенетического дерева с помощью iqtree

from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus
from ete3 import Tree
import subprocess

from config import MSA_MATRIX_PATH, TREES_NUM, TREE_PATH, OUTGROUP
from support_functions import mkdir

# Проверка и анализ топологии локальных деревьев
def check_topology(trees):
    topologies_found = {}

    for tree in trees:
        is_unique = True

        for topology in topologies_found:
            compare = tree.compare(topology, unrooted=True)
            rf = compare["rf"]

            if rf == 0:
                topologies_found[topology] += 1
                is_unique = False
        
        if is_unique:
            topologies_found[tree] = 1
    
    print(f"Different topologies among local trees: {len(topologies_found)}")
    print("Details:")

    for topology in topologies_found:
        percentage = topologies_found[topology] * 100 / TREES_NUM

        print(f"{topology} - {topologies_found[topology]} trees ({percentage}%)")

# Построение консенсусного дерева
def build_majority_tree(trees):
    majority_tree = majority_consensus(trees, cutoff=0.5)

    majority_tree_file = f"{TREE_PATH}/majority_tree.nwk"
    Phylo.write(majority_tree, majority_tree_file, "newick")

    print(f"Consensus majority tree was saved in {majority_tree_file}, it looks like:")
    Phylo.draw_ascii(majority_tree)

# Запуск IQ-TREE
def run_iqtree():
    mkdir(TREE_PATH)

    trees_phylo = []
    trees_ete3 = []

    for cnt in range(TREES_NUM):
        matrix_fasta = f"{MSA_MATRIX_PATH}/matrix_{cnt}.fasta"

        cmd = f"iqtree2 -s {matrix_fasta} -o {OUTGROUP} -pre {TREE_PATH}/Tree_{cnt} -B 1000 --scfl 100 -m MFP+ASC -T 8 -mem 8G"
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        treefile = f"{TREE_PATH}/Tree_{cnt}.treefile"

        trees_phylo.append(Phylo.read(treefile, "newick"))
        trees_ete3.append(Tree(treefile))
    
    build_majority_tree(trees_phylo)
    check_topology(trees_ete3)