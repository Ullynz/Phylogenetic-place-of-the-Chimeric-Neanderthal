# Строительство филогенетического дерева с помощью iqtree

from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus
import subprocess

from config import MSA_MATRIX_PATH, TREES_NUM, TREE_PATH, OUTGROUP
from support_functions import mkdir

def run_iqtree():
    mkdir(TREE_PATH)

    trees = []

    for cnt in range(TREES_NUM):
        matrix_fasta = f"{MSA_MATRIX_PATH}/matrix_{cnt}.fasta"

        cmd = f"iqtree2 -s {matrix_fasta} -o {OUTGROUP} -pre {TREE_PATH}/Tree_{cnt} -B 1000 --scfl 100 -m MFP+ASC -T 8 -mem 8G"
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        treefile = f"{TREE_PATH}/Tree_{cnt}.treefile"
        trees.append(Phylo.read(treefile, "newick"))

    majority_tree = majority_consensus(trees, cutoff=0.5)

    majority_tree_file = f"{TREE_PATH}/majority_tree.nwk"
    Phylo.write(majority_tree, majority_tree_file, "newick")

    print(f"Consensus majority tree was saved in {majority_tree_file}, it looks like:")
    Phylo.draw_ascii(majority_tree)