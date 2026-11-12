#Helpers for FASTQ filtering (streaming, stdlib only).

from typing import Iterator, Tuple


def parse_bounds(bounds, default_low: float = 0.0) -> Tuple[float, float]:
    """
    Превращает одиночное число в кортеж (0, число), а кортеж оставляет как есть.
    """
    if isinstance(bounds, (int, float)):
        return (float(default_low), float(bounds))
    low, high = bounds
    return (float(low), float(high))


def read_fastq(input_fastq: str) -> Iterator[Tuple[str, str, str, str]]:
    """
    Генератор (ленивое чтение). Выдает по одному риду за раз (4 строки).
    Благодаря этому мы не загружаем весь файл в память.
    """
    with open(input_fastq, "r") as f:
        while True:
            header = f.readline()
            if not header:  # Конец файла
                return
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            # Проверка на корректность формата
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record at header: {header!r}")
                
            yield header.rstrip("\n"), seq, plus, qual


def is_gc_in_bounds(seq: str, low: float, high: float) -> bool:
    """Проверяет, что GC-состав в процентах входит в [low, high]."""
    if not seq:
        return False
    # Считаем G и C (в любом регистре), делим на длину и умножаем на 100
    gc = sum(1 for c in seq if c in "GCgc") / len(seq) * 100
    return low <= gc <= high


def is_length_in_bounds(seq: str, low: float, high: float) -> bool:
    """Проверяет, что длина последовательности входит в [low, high]."""
    return low <= len(seq) <= high


def is_quality_ok(qual: str, threshold: float) -> bool:
    """Проверяет, что среднее качество Phred33 >= threshold"""
    if not qual:
        return False
    # ord(c) - 33 переводит ASCII-символ в числовое значение качества
    mean_q = sum(ord(c) - 33 for c in qual) / len(qual)
    return mean_q >= threshold


def passes_filters(
    seq: str,
    qual: str,
    gc_bounds: Tuple[float, float],
    length_bounds: Tuple[float, float],
    quality_threshold: float,
) -> bool:
    """
    Объединяет все три проверки в одну. 
    Возвращает True, только если рид прошел ВСЕ фильтры.
    """
    return (
        is_gc_in_bounds(seq, *gc_bounds)
        and is_length_in_bounds(seq, *length_bounds)
        and is_quality_ok(qual, quality_threshold)
    )