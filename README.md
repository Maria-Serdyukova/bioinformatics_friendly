# bioinformatics_friendly

Биоинформатические утилиты: фильтрация FASTQ, конвертация FASTA и парсинг BLAST.

## Установка
Сторонних зависимостей нет. Требуется Python 3.10+.

## Использование

```python
from main import filter_fastq, convert_multiline_fasta_to_oneline, parse_blast_output

# FASTQ фильтрация (сохраняет в папку filtered/)
filter_fastq("input.fastq", "output.fastq", gc_bounds=(20, 80), quality_threshold=30)

# FASTA (многострочный -> однострочный)
convert_multiline_fasta_to_oneline("input.fasta")  # создаст input_oneline.fasta

# BLAST (извлечь лучшие совпадения)
parse_blast_output("blast.txt", "proteins.txt")