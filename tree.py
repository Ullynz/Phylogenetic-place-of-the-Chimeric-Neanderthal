# Строительство филогенетического дерева с помощью iqtree

import subprocess

from config import PROJ_PATH
from support_functions import check_file, mkdir

def run_iqtree(MSA_matrix):
    tree_dir = f"{PROJ_PATH}/trees"
    mkdir(tree_dir)

    # В файле num.txt храню индекс последнего построенного дерева, потом использую это число в префиксе для выходных файлов iqtree
    tree_num_file = f"{tree_dir}/num.txt"
    if not check_file(tree_num_file):
        with open(tree_num_file, "w") as f:
            f.write("0")
    
    with open(tree_num_file, 'r') as f:
        nxt_num = int(f.read()) + 1
    
    with open(tree_num_file, "w") as f:
        f.write(str(nxt_num))

    cmd = f"iqtree2 -s {MSA_matrix} -pre {tree_dir}/FinalTree_{nxt_num} -B 1000 -T 8 -mem 8G"
    subprocess.run(cmd, shell=True)