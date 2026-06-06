# Пути к файлам, необходимым для пайплайна:

# Директории с результатами работы DAIseg
DAISEG_PATHS = [
    "/home/usnasonova/DAIseg/TSI.YRI",
    "/home/usnasonova/DAIseg/IBS.YRI",
    "/home/usnasonova/DAIseg/CHB"
]

# json-ы для DAIseg
DAISEG_JSONS = [
    "/home/usnasonova/DAIseg/TSI.YRI/jsons/TSI.YRI.grch37.chr21.json",
    "/home/usnasonova/DAIseg/IBS.YRI/jsons/IBS.YRI.grch37.chr21.json",
    ""
]

# Пути к output-ам DAIseg
def get_daiseg_output(daiseg_path, population, chrom):
    return f"{daiseg_path}/{population}.YRI.grch37.chr{chrom}.em.tsv"

# Пути к coverage
def get_modern_coverage(chrom):
    return f"/home/share/human.data/1000GP/1000GP.grch37/bed/chr{chrom}.renamed.bed"

def get_altai_coverage(chrom):
    return f"/home/share/human.data/neand/altai/bed/chr{chrom}_mask.bed.gz"

def get_vindija_coverage(chrom):
    return f"/home/share/human.data/neand/33.19/bed/chr{chrom}_mask.bed.gz"

def get_outgroup_coverage(chrom):
    return f"/home/usnasonova/project/Ust/masks/chr{chrom}.Map35_99.MQ30.Cov.indels.TRF.bed.gz"
    #return get_modern_coverage(chrom)

# Пути к vcf 
def get_modern_pop_vcf(chrom):
    return f"/home/share/human.data/1000GP/1000GP.grch37/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

def get_vindija_vcf(chrom):
    return f"/home/share/human.data/neand/33.19/chr{chrom}_mq25_mapab100.vcf.gz"

def get_altai_vcf(chrom):
    return f"/home/share/human.data/neand/altai/chr{chrom}_mq25_mapab100.vcf.gz"

def get_outgroup_vcf(chrom):
    return f"/home/usnasonova/project/Ust/Ust_Ishim.hg19_1000g.{chrom}.mod.vcf.gz"
    #return get_modern_pop_vcf(chrom)


# Параметры самого пайплайна:

# Рабочая директория
PROJ_PATH = "/home/usnasonova/project" 

# Популяции
MODERN_POP = ["TSI", "IBS", "CHB"]
NEAND_POP = ["Vindija", "Altai"]
OUTGROUP = "Ust"

# Директория, в которой появятся маски для пересеченных архаичных сегментов
ARCHAIC_SECTIONS_PATH = f"{PROJ_PATH}/intersected_archaic_sections"

# Директория для fasta файлов MSA матриц
MSA_MATRIX_PATH = f"{PROJ_PATH}/matrix"

# Директория для результатов iqtree
TREE_PATH = f"{PROJ_PATH}/trees"

# Число локальных деревьев
TREES_NUM = 100

# Размер окон, по которым строятся деревья (-1, если нужен максимальный, чтобы TREE_NUM деревьев покрывали все позиции)
WINDOW_SIZE = -1

