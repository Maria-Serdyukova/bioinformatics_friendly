#Helpers for reading bioinformatics files (FASTA, BLAST)

import os
from helpers.fastq import parse_bounds, read_fastq, passes_filters
from helpers.bio_files import iter_fasta_records, default_fasta_name, extract_first_descriptions


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2 ** 32),
    quality_threshold: float = 0,
) -> None:
    
    """
    Фильтрует FASTQ-риды по GC-составу, длине и среднему качеству Phred33.
    Результат сохраняется в папку 'filtered/'.
    
    Аргументы:
    input_fastq: путь до входного FASTQ файла
    output_fastq: имя выходного файла (сохраняется в 'filtered/')
    gc_bounds: (low, high) в процентах или одно число (верхняя граница)
    length_bounds: (low, high) или одно число (верхняя граница)
    quality_threshold: минимальное среднее качество (включительно)
    """
    # 1. Нормализуем границы (если передано одно число, сделает (0, число))
    gc = parse_bounds(gc_bounds, 0)
    length = parse_bounds(length_bounds, 0)

    # 2. Создаем папку filtered и формируем путь
    os.makedirs("filtered", exist_ok=True)
    out_path = os.path.join("filtered", output_fastq)

    # 3. Открываем файл и читаем по одному риду
    with open(out_path, "w") as out:
        for header, seq, plus, qual in read_fastq(input_fastq):
            # Если рид прошел все проверки — записываем
            if passes_filters(seq, qual, gc, length, quality_threshold):
                out.write(f"{header}\n{seq}\n{plus}\n{qual}\n")


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str | None = None
) -> None:
    
    """
    Конвертирует многострочный FASTA в FASTA с одной строкой на последовательность.
    
    Аргументы:
    input_fasta: путь до входного FASTA
    output_fasta: путь до выходного (если None, создастся '<input>_oneline.fasta')
    """
    
    # Если имя не передали — генерируем его
    if output_fasta is None:
        output_fasta = default_fasta_name(input_fasta)
        
    with open(output_fasta, "w") as out:
        # iter_fasta_records склеивает куски последовательности в одну строку
        for header, seq in iter_fasta_records(input_fasta):
            out.write(f"{header}\n{seq}\n")


def parse_blast_output(input_file: str, output_file: str) -> None:
    
    """
    Достает первое совпадение (Description) для каждого QUERY из BLAST txt.
    Сохраняет отсортированный список в один столбец.
    """
    
    # Получаем отсортированный список названий
    descriptions = extract_first_descriptions(input_file)
    
    # Записываем в файл
    with open(output_file, "w") as out:
        for d in descriptions:
            out.write(d + "\n")