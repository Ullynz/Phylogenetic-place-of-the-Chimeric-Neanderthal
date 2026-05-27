# Строительство филогенетического дерева с помощью iqtree

import subprocess

from config import MSA_MATRIX_PATH, TREES_NUM, TREE_PATH
from support_functions import mkdir

def run_iqtree():
    mkdir(TREE_PATH)

    for cnt in range(TREES_NUM):
        matrix_fasta = f"{MSA_MATRIX_PATH}/matrix_{cnt}.fasta"

        cmd = f"iqtree2 -s {matrix_fasta} -pre {TREE_PATH}/Tree_{cnt} -B 1000 --scfl 100 -m MFP+ASC -T 8 -mem 8G"
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
