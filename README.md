# algorithms_hw1

# Инструкция по использованию программы

## Основные аргументы командной строки

| Аргумент | Описание | Тип | Обязательный |
|----------|----------|-----|--------------|
| `--genome1` | FASTA-файл первого генома (низкий GC) | путь к файлу | Да |
| `--genome2` | FASTA-файл второго генома (высокий GC) | путь к файлу | Да |
| `--chimera` | Флаг генерации химерной последовательности | флаг | Нет |
| `--length` | Длина химерной последовательности (по умолчанию: 10000) | число | Нет |
| `--avg-fragment` | Средняя длина фрагментов (по умолчанию: 300) | число | Нет |
| `--input` | Входная последовательность для декодирования | путь к файлу | Нет |
| `--true-labels` | Файл с истинными метками для оценки точности | путь к файлу | Нет |
| `--output` | Файл для сохранения предсказаний | путь к файлу | Нет |
| `--verbose` | Подробный вывод информации | флаг | Нет |

## Установка зависимостей

```bash
pip install biopython
```

## Скачивание геномов
необходимо скачать геномы, я это делала вручную с NCBI, но можно через ```wget```

## Генерация химерной последовательности

```bash
python3 hw1_algorithms_Kosareva.py \
    --genome1 Rickettsia_prowazekii.fasta \
    --genome2 Deinococcus_radiodurans.fasta \
    --chimera --length 10000 --avg-fragment 300 --verbose
```

## Анализ химерной последовательности

```bash
python3 hw1_algorithms_Kosareva.py \
    --genome1 Rickettsia_prowazekii.fasta \
    --genome2 Deinococcus_radiodurans.fasta \
    --input chimera.fasta \
    --true-labels chimera_true_labels.txt \
    --output chimera_predictions.txt --verbose
```

## Анализ чистых геномов

```bash
python3 hw1_algorithms_Kosareva.py \
    --genome1 Rickettsia_prowazekii.fasta \
    --genome2 Deinococcus_radiodurans.fasta \
    --input test_genome1.fasta \
    --true-labels test_genome1_true_labels.txt \
    --output genome1_predictions.txt

python3 hw1_algorithms_Kosareva.py \
    --genome1 Rickettsia_prowazekii.fasta \
    --genome2 Deinococcus_radiodurans.fasta \
    --input test_genome2.fasta \
    --true-labels test_genome2_true_labels.txt \
    --output genome2_predictions.txt
```
