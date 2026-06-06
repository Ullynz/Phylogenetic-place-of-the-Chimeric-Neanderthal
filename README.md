# Филогенетическое место «химерного» неандертальца

## Структура проекта

`main.py` - главная функция

`masks.py` - формирование масок с участками неандертальской интрогрессии и хорошим покрытием

`filter_vcf.py` - фильтрация vcf по найденным маскам

`chimeric_genome.py` - построение "химерного" генома

`msa_matrix.py` - построение MSA матриц по химерным геномам

`tree.py` - построение филогенетических деревьев, формирование консенсусного дерева

`support_functions.py` - вспомогательные функции

`config.py` - файл с конфигурациями

### Конфигурации в `config.py`

Пути, связанные с DAIseg:
- `DAISEG_PATHS`: пути к директориям с результатами работы DAIseg для современных популяций
- `DAISEG_JSONS`: config json-ы, поданные в DAIseg

Пути для текущего пайплайна:
- `PROJ_PATH`: рабочая директория для выполнения пайплайна
- `ARCHAIC_SECTIONS_PATH`: выходная директория, в которой появятся маски для пересеченных архаичных сегментов
- `MSA_MATRIX_PATH`: выходная директория для fasta файлов MSA матриц
- `TREE_PATH`: выходная директория для результатов iqtree

Параметры пайплайна:
- `MODERN_POP`: список современных популяций
- `NEAND_POP = ["Vindija", "Altai"]`: список неандертальских популяций
- `OUTGROUP`: внешняя группа
- `TREES_NUM`: число строимых локальных деревьев
- `WINDOW_SIZE`: размер окон, по которым строятся деревья (-1, если нужен максимальный, чтобы TREE_NUM деревьев покрывали все позиции)

Функции для нахождения путей к различным данным по хромосоме:
- `get_daiseg_output(daiseg_path, population, chrom)`: путь к output-у DAIseg для популяции `population` хромосомы `chrom` и заданного пути к директории с результатами
- `get_modern_coverage(chrom)`: путь к coverage маске для современных популяций
- `get_altai_coverage(chrom)`: путь к coverage маске для Altai
- `get_vindija_coverage(chrom)`: путь к coverage маске для Vindija
- `get_outgroup_coverage(chrom)`: путь к coverage маске для outgroup
- `get_modern_pop_vcf(chrom)`: путь к vcf современных популяций
- `get_vindija_vcf(chrom)`: геном Vindija
- `get_altai_vcf(chrom)`: геном Altai
- `get_outgroup_vcf(chrom)`: геном outgroup-ы

## Использование
Для запуска пайплайна сперва необходимо для каждой из рассматриваемых современных популяций запустить [`DAIseg`](https://github.com/Genomics-HSE/DAIseg.git). Также нам потребуются геномные данные (vcf и coverage-маски) для outgroup-ы - если ранее они не были использованы для `DAIseg`, их необходимо предварительно скачать.

Далее необходимо склонировать данный репозиторий, сделать это можно следующим образом:
```
git clone https://github.com/Ullynz/Phylogenetic-place-of-the-Chimeric-Neanderthal.git
```
Задать необходимую конфигурацию в `config.py`, и затем запустить пайплайн как:
```
python3 main.py
```
