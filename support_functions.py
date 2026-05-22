# Модуль с вспомогательными функциями

from pathlib import Path
import random

def rnd():
    return int(random.random() >= 0.5)

def check_file(file_name):
    file_path = Path(file_name)

    if file_path.exists() and file_path.stat().st_size > 0:
        return True
    return False

def mkdir(dir_name):
    Path(dir_name).mkdir(parents=True, exist_ok=True)