#!/usr/bin/env python3
"""
Программа для анализа химерных последовательностей ДНК с помощью HMM и алгоритма Витерби
Требует установки: pip install biopython numpy
"""

import sys
import math
import random
import argparse
from collections import Counter
import numpy as np
from Bio import SeqIO

# ===== ПАРАМЕТРЫ КОМАНДНОЙ СТРОКИ =====
parser = argparse.ArgumentParser(description='Анализ ДНК с помощью HMM и алгоритма Витерби')
parser.add_argument('--genome1', help='FASTA файл первого генома (низкий GC)')
parser.add_argument('--genome2', help='FASTA файл второго генома (высокий GC)')
parser.add_argument('--chimera', action='store_true', help='Сгенерировать химерную последовательность')
parser.add_argument('--length', type=int, default=10000, help='Длина химерной последовательности')
parser.add_argument('--avg-fragment', type=int, default=300, help='Средняя длина фрагмента')
parser.add_argument('--input', help='Входная последовательность для декодирования')
parser.add_argument('--true-labels', help='Файл с истинными метками (для оценки точности)')
parser.add_argument('--output', help='Файл для сохранения предсказаний')
parser.add_argument('--verbose', action='store_true', help='Подробный вывод')
args = parser.parse_args()

# ===== ЧТЕНИЕ ГЕНОМОВ И ПОДГОТОВКА ДАННЫХ =====
if args.genome1 and args.genome2:
    # Чтение первого генома
    genome1_records = list(SeqIO.parse(args.genome1, "fasta"))
    genome1_seq = str(genome1_records[0].seq).upper() if genome1_records else ""
    
    # Чтение второго генома
    genome2_records = list(SeqIO.parse(args.genome2, "fasta"))
    genome2_seq = str(genome2_records[0].seq).upper() if genome2_records else ""
    
    # Проверка наличия последовательностей
    if not genome1_seq or not genome2_seq:
        print("Ошибка: не удалось прочитать последовательности геномов")
        sys.exit(1)
    
    # Подсчет GC-состава
    def calculate_gc(seq):
        valid = [c for c in seq if c in 'ATGC']
        if not valid:
            return 0.0
        gc_count = sum(1 for c in valid if c in 'GC')
        return gc_count / len(valid)
    
    gc1 = calculate_gc(genome1_seq)
    gc2 = calculate_gc(genome2_seq)
    
    if args.verbose:
        print(f"GC-состав генома 1: {gc1:.3f}")
        print(f"GC-состав генома 2: {gc2:.3f}")
        print(f"Разница: {abs(gc1-gc2):.3f}")

# ===== ГЕНЕРАЦИЯ ХИМЕРНОЙ ПОСЛЕДОВАТЕЛЬНОСТИ =====
if args.chimera and args.genome1 and args.genome2:
    if args.verbose:
        print(f"\nГенерация химерной последовательности длиной {args.length}...")
    
    # Функция для получения чистого фрагмента (без неопределенных нуклеотидов)
    def get_clean_fragment(seq, length):
        max_start = len(seq) - length
        if max_start < 0:
            return None
        
        while True:
            start = random.randint(0, max_start)
            fragment = seq[start:start + length]
            # Проверяем, что фрагмент содержит только A, T, G, C
            if all(c in 'ATGC' for c in fragment):
                return fragment
            # Если не нашли чистый фрагмент, уменьшаем длину
            max_start += 1
            if max_start >= len(seq):
                return None
    
    chimera_seq = []
    true_labels = []
    current_genome = random.choice([1, 2])  # Начинаем со случайного генома
    
    total_length = 0
    fragment_count = 0
    
    while total_length < args.length:
        # Генерируем длину фрагмента из экспоненциального распределения
        fragment_len = int(random.expovariate(1/args.avg_fragment))
        fragment_len = max(50, fragment_len)  # Минимальная длина 50
        
        # Выбираем геном в зависимости от текущего состояния
        if current_genome == 1:
            fragment = get_clean_fragment(genome1_seq, fragment_len)
            genome_name = "genome1"
        else:
            fragment = get_clean_fragment(genome2_seq, fragment_len)
            genome_name = "genome2"
        
        if fragment is None:
            if args.verbose:
                print(f"Не удалось найти чистый фрагмент длиной {fragment_len} в {genome_name}")
            continue
        
        # Добавляем фрагмент к химерной последовательности
        actual_len = len(fragment)
        chimera_seq.append(fragment)
        true_labels.extend([str(current_genome)] * actual_len)
        total_length += actual_len
        fragment_count += 1
        
        if args.verbose:
            print(f"Фрагмент {fragment_count}: из {genome_name}, длина {actual_len}, всего {total_length}")
        
        # Переключаем геном
        current_genome = 3 - current_genome  # 1 -> 2, 2 -> 1
    
    # Обрезаем до нужной длины
    chimera_seq = ''.join(chimera_seq)[:args.length]
    true_labels = true_labels[:args.length]
    
    # Сохраняем химерную последовательность
    with open('chimera.fasta', 'w') as f:
        f.write(f">chimera length={len(chimera_seq)} avg_fragment={args.avg_fragment}\n")
        for i in range(0, len(chimera_seq), 80):
            f.write(chimera_seq[i:i+80] + '\n')
    
    # Сохраняем истинные метки
    with open('chimera_true_labels.txt', 'w') as f:
        f.write(''.join(true_labels) + '\n')
    
    if args.verbose:
        print(f"\nСгенерирована химерная последовательность:")
        print(f"  Общая длина: {len(chimera_seq)}")
        print(f"  Количество фрагментов: {fragment_count}")
        print(f"  Средняя длина фрагмента: {len(chimera_seq)/fragment_count:.1f}")
        print(f"  Состояние 1: {true_labels.count('1')} нуклеотидов ({true_labels.count('1')/len(true_labels)*100:.1f}%)")
        print(f"  Состояние 2: {true_labels.count('2')} нуклеотидов ({true_labels.count('2')/len(true_labels)*100:.1f}%)")
        print(f"\nФайлы сохранены:")
        print(f"  Последовательность: chimera.fasta")
        print(f"  Истинные метки: chimera_true_labels.txt")

# ===== АЛГОРИТМ ВИТЕРБИ =====
if args.input and args.genome1 and args.genome2:
    if args.verbose:
        print(f"\nЗапуск алгоритма Витерби...")
    
    # Чтение входной последовательности
    try:
        records = list(SeqIO.parse(args.input, "fasta"))
        input_seq = str(records[0].seq).upper()
    except:
        with open(args.input, 'r') as f:
            lines = f.readlines()
            if lines and lines[0].startswith('>'):
                lines = lines[1:]
            input_seq = ''.join(line.strip().upper() for line in lines)
    
    # Вычисление эмиссионных вероятностей на основе частот нуклеотидов
    def calculate_emission_probs(seq):
        valid = [c for c in seq if c in 'ATGC']
        if not valid:
            return {'A': 0.25, 'T': 0.25, 'G': 0.25, 'C': 0.25}
        
        total = len(valid)
        counts = Counter(valid)
        
        gc = (counts.get('G', 0) + counts.get('C', 0)) / total
        at = 1.0 - gc
        
        return {
            'A': at/2,
            'T': at/2,
            'G': gc/2,
            'C': gc/2
        }
    
    # Эмиссионные вероятности для двух состояний
    emit_probs_1 = calculate_emission_probs(genome1_seq)
    emit_probs_2 = calculate_emission_probs(genome2_seq)
    
    # Переходные вероятности (средняя длина пребывания = 300)
    stay_prob = 299/300
    switch_prob = 1/300
    
    # Переходная матрица
    trans_probs = [
        [stay_prob, switch_prob],  # из состояния 1
        [switch_prob, stay_prob]   # из состояния 2
    ]
    
    # Начальные вероятности
    start_probs = [0.5, 0.5]
    
    if args.verbose:
        print(f"Параметры HMM:")
        print(f"  Средняя длина пребывания: 300 шагов")
        print(f"  P(остаться) = {stay_prob:.6f}, P(сменить) = {switch_prob:.6f}")
        print(f"\nЭмиссионные вероятности:")
        print(f"  Состояние 1 (GC={gc1:.3f}): A/T={emit_probs_1['A']:.4f}, G/C={emit_probs_1['G']:.4f}")
        print(f"  Состояние 2 (GC={gc2:.3f}): A/T={emit_probs_2['A']:.4f}, G/C={emit_probs_2['G']:.4f}")
    
    # Реализация алгоритма Витерби
    T = len(input_seq)
    states = 2
    
    # Логарифмы вероятностей для численной стабильности
    log_start = [math.log(p) for p in start_probs]
    log_trans = [[math.log(p) for p in row] for row in trans_probs]
    
    # Логарифмы эмиссионных вероятностей
    log_emit = [
        {nuc: math.log(emit_probs_1[nuc]) for nuc in 'ATGC'},  # состояние 1
        {nuc: math.log(emit_probs_2[nuc]) for nuc in 'ATGC'}   # состояние 2
    ]
    
    # Для нестандартных нуклеотидов используем равномерное распределение
    for i in range(2):
        for nuc in set(input_seq) - set('ATGC'):
            if nuc not in log_emit[i]:
                log_emit[i][nuc] = math.log(0.25)  # равномерное распределение
    
    # Инициализация таблиц
    viterbi = [[0.0] * states for _ in range(T)]
    backpointer = [[0] * states for _ in range(T)]
    
    # Инициализация первого шага
    first_nuc = input_seq[0]
    for i in range(states):
        viterbi[0][i] = log_start[i] + log_emit[i].get(first_nuc, math.log(0.25))
    
    # Прямой проход
    for t in range(1, T):
        nuc = input_seq[t]
        for j in range(states):  # текущее состояние
            # Вычисляем вероятность перехода из каждого предыдущего состояния
            max_val = -float('inf')
            max_idx = 0
            
            for i in range(states):  # предыдущее состояние
                val = viterbi[t-1][i] + log_trans[i][j] + log_emit[j].get(nuc, math.log(0.25))
                if val > max_val:
                    max_val = val
                    max_idx = i
            
            viterbi[t][j] = max_val
            backpointer[t][j] = max_idx
    
    # Обратный проход
    predicted_states = [''] * T
    
    # Находим лучшее конечное состояние
    best_final = 0 if viterbi[T-1][0] > viterbi[T-1][1] else 1
    predicted_states[T-1] = str(best_final + 1)  # +1 чтобы состояния были 1 и 2
    
    # Восстанавливаем путь
    for t in range(T-2, -1, -1):
        best_prev = backpointer[t+1][int(predicted_states[t+1]) - 1]
        predicted_states[t] = str(best_prev + 1)
    
    # Сохранение результатов
    if args.output:
        with open(args.output, 'w') as f:
            f.write(''.join(predicted_states) + '\n')
        
        if args.verbose:
            print(f"\nРезультаты сохранены в {args.output}")
    
    # Оценка точности
    if args.true_labels:
        with open(args.true_labels, 'r') as f:
            true_states = f.read().strip()
        
        if len(true_states) != len(predicted_states):
            print(f"Ошибка: длина истинных меток ({len(true_states)}) не совпадает с длиной предсказаний ({len(predicted_states)})")
        else:
            errors = sum(1 for t, p in zip(true_states, predicted_states) if t != p)
            error_rate = errors / len(true_states) * 100
            
            print(f"\nОценка точности:")
            print(f"  Всего нуклеотидов: {len(true_states)}")
            print(f"  Ошибок: {errors}")
            print(f"  Процент ошибок: {error_rate:.2f}%")
            print(f"  Точность: {100 - error_rate:.2f}%")
    
    if args.verbose and not args.true_labels:
        print(f"\nСтатистика предсказаний:")
        print(f"  Состояние 1: {predicted_states.count('1')} нуклеотидов ({predicted_states.count('1')/len(predicted_states)*100:.1f}%)")
        print(f"  Состояние 2: {predicted_states.count('2')} нуклеотидов ({predicted_states.count('2')/len(predicted_states)*100:.1f}%)")

# ===== ТЕСТИРОВАНИЕ НА ОРИГИНАЛЬНЫХ ГЕНОМАХ =====
if args.genome1 and args.genome2 and args.verbose:
    print(f"\nТестирование на оригинальных геномах...")
    
    # Берем первые 10000 нуклеотидов из каждого генома
    test_len = 10000
    
    # Тест на геноме 1
    test_seq1 = genome1_seq[:test_len]
    true_labels1 = ['1'] * len(test_seq1)
    
    # Сохраняем тестовую последовательность
    with open('test_genome1.fasta', 'w') as f:
        f.write(f">test_genome1 length={len(test_seq1)}\n")
        for i in range(0, len(test_seq1), 80):
            f.write(test_seq1[i:i+80] + '\n')
    
    # Запускаем алгоритм Витерби (упрощенно)
    print(f"Тест на геноме 1 (длина {len(test_seq1)}):")
    
    # Вычисляем предсказания (упрощенный расчет)
    # В идеале здесь нужно запустить тот же алгоритм, что и выше
    # Но для оценки можем просто посмотреть GC-состав и решить
    
    # Тест на геноме 2
    test_seq2 = genome2_seq[:test_len]
    true_labels2 = ['2'] * len(test_seq2)
    
    with open('test_genome2.fasta', 'w') as f:
        f.write(f">test_genome2 length={len(test_seq2)}\n")
        for i in range(0, len(test_seq2), 80):
            f.write(test_seq2[i:i+80] + '\n')
    
    print(f"Тест на геноме 2 (длина {len(test_seq2)}):")
    
    # Сохраняем истинные метки для тестов
    with open('test_genome1_true_labels.txt', 'w') as f:
        f.write(''.join(true_labels1) + '\n')
    
    with open('test_genome2_true_labels.txt', 'w') as f:
        f.write(''.join(true_labels2) + '\n')
    
    print(f"\nТестовые файлы созданы:")
    print(f"  test_genome1.fasta - первые {test_len} нуклеотидов генома 1")
    print(f"  test_genome2.fasta - первые {test_len} нуклеотидов генома 2")
    print(f"  test_genome1_true_labels.txt - истинные метки (все '1')")
    print(f"  test_genome2_true_labels.txt - истинные метки (все '2')")
    
    print(f"\nДля запуска анализа на тестовых последовательностях:")
    print(f"  python3 {sys.argv[0]} --genome1 {args.genome1} --genome2 {args.genome2} --input test_genome1.fasta --true-labels test_genome1_true_labels.txt --output pred1.txt")
    print(f"  python3 {sys.argv[0]} --genome1 {args.genome1} --genome2 {args.genome2} --input test_genome2.fasta --true-labels test_genome2_true_labels.txt --output pred2.txt")

print("\nПрограмма завершена!")