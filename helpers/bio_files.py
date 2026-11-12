#Helpers for reading bioinformatics files (FASTA, BLAST).

import re
from typing import Iterator, List, Tuple, Union


def iter_fasta_records(input_fasta: str) -> Iterator[Tuple[str, str]]:
    """
    Лениво читает FASTA и выдает (заголовок, последовательность).
    Склеивает многострочные последовательности в одну строку.
    """
    header: Union[str, None] = None
    chunks: List[str] = []
    
    with open(input_fasta, "r") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                # Если нашли новый заголовок, отдаем предыдущий
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line, []
            elif line:
                chunks.append(line.strip())
                
    # Отдаем последнюю последовательность в конце файла
    if header is not None:
        yield header, "".join(chunks)


def default_fasta_name(input_fasta: str) -> str:
    """ Генерирует имя выходного файла: 'input.fasta' -> 'input_oneline.fasta' """
    base = input_fasta.replace("\\", "/").rsplit("/", 1)[-1]
    # Убираем расширение .fa, .fasta, .fna
    stem = re.sub(r"\.(fa|fasta|fna)$", "", base, flags=re.IGNORECASE)
    return f"{stem}_oneline.fasta"


def extract_first_descriptions(blast_txt: str) -> List[str]:
    """
    Парсит BLAST txt. Ищет заголовки 'Sequences producing...'
    Берет первую строку данных (лучшее совпадение) и вытаскивает Description.
    """
    results: List[str] = []
    in_block = False
    header_seen = False
    
    with open(blast_txt, "r") as f:
        for line in f:
            stripped = line.strip()
            
            # Вход в новый блок результатов
            if "Sequences producing significant alignments" in stripped:
                in_block = True
                header_seen = False
                continue
                
            if not in_block or not stripped:
                continue
                
            # Пропускаем шапку таблицы 
            if not header_seen:
                header_seen = True
                continue
                
            # Это первая строка с данными, разбиваем по 2+ пробелам
            parts = re.split(r"\s{2,}", stripped)
            if parts:
                results.append(parts[0].strip())  # Берем столбец Description
                
            # Выключаем блок, так как нам нужна только первая строка
            in_block = False
            
    # Сортируем и убираем дубликаты (set)
    return sorted(set(results))

    #+ в этой папке есть helpers/__init__.py - нужен для того, чтобы Python понимал, что helpers — это пакет.