# Пути к файлам, необходимым для пайплайна:

# Директории с результатами работы DAIseg
DAISEG_POP_A_PATH = "/home/usnasonova/DAIseg/TSI.YRI"
DAISEG_POP_B_PATH = "/home/usnasonova/DAIseg/IBS.YRI"

# json-ы для DAIseg
DAISEG_POP_A_JSON = "/home/usnasonova/DAIseg/TSI.YRI/jsons/TSI.YRI.grch37.chr21.json"
DAISEG_POP_B_JSON = "/home/usnasonova/DAIseg/IBS.YRI/jsons/IBS.YRI.grch37.chr21.json"

# Пути к output-ам DAIseg
def get_daiseg_pop_A_output(chrom):
    return f"{DAISEG_POP_A_PATH}/TSI.YRI.grch37.chr{chrom}.em.tsv"

def get_daiseg_pop_B_output(chrom):
    return f"{DAISEG_POP_B_PATH}/IBS.YRI.grch37.chr{chrom}.em.tsv"

# Пути к coverage
def get_daiseg_coverage(daiseg_pop_path, chrom):
    return f"{daiseg_pop_path}/coverage_1kG.chr{chrom}.grch37.bed"

def get_altai_coverage(chrom):
    return f"/home/share/human.data/neand/altai/bed/chr{chrom}_mask.bed.gz"

def get_vindija_coverage(chrom):
    return f"/home/share/human.data/neand/33.19/bed/chr{chrom}_mask.bed.gz"

def get_outgroup_coverage(chrom):
    return f"/home/usnasonova/project/Ust/masks/chr{chrom}.Map35_99.MQ30.Cov.indels.TRF.bed.gz"

# Путь к фильтрованным после DAIseg vcf 
def get_modern_pop_vcf(daiseg_pop_path, chrom):
    return f"{daiseg_pop_path}/1kG_filtered.chr{chrom}.grch37.vcf.gz"

# Геном Vindija
def get_vindija_vcf(chrom):
    return f"/home/share/human.data/neand/33.19/chr{chrom}_mq25_mapab100.vcf.gz"

# Геном Altai
def get_altai_vcf(chrom):
    return f"/home/share/human.data/neand/altai/chr{chrom}_mq25_mapab100.vcf.gz"

# Геном outgroup-ы
def get_outgroup_vcf(chrom):
    return f"/home/usnasonova/project/Ust/Ust_Ishim.hg19_1000g.{chrom}.mod.vcf.gz"


# Параметры самого пайплайна:

# Рабочая директория
PROJ_PATH = "/home/usnasonova/project" 

# Популяции
MODERN_POP = ["TSI", "IBS"]
NEAND_POP = ["Vindija", "Altai"]
OUTGROUP = "Ust"

# Директория, в которой появятся маски для пересеченных архаичных сегментов
ARCHAIC_SECTIONS_PATH = f"{PROJ_PATH}/intersected_archaic_sections"

# Директория для fasta файлов MSA матриц
MSA_MATRIX_PATH = f"{PROJ_PATH}/matrix"

# Директория для результатов iqtree
TREE_PATH = f"{PROJ_PATH}/trees"

# Число локальных деревьев
TREES_NUM = 40

# Размер окон, по которым строятся деревья (-1, если нужен максимальный, чтобы TREE_NUM деревьев покрывали все позиции)
WINDOW_SIZE = -1

