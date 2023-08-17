import argparse
import os
import re
import sys
from collections import Counter, defaultdict

import matplotlib.pyplot as plt
import numpy as np


def progress(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s  %s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()  # As suggested by Rom Ruben



argParser = argparse.ArgumentParser()
argParser.add_argument("-f", "--filenames", type=str, required=True, help="Input .sam files. Example: miRNA_S7872Nr2.1.fastq.gz.sam,miRNA_S7872Nr3.1.fastq.gz,miRNA_S7872Nr4.1.fastq.gz")
argParser.add_argument("-c", "--count", type=int, default=1, help="Сount nucleotides in start of a read")
argParser.add_argument("-o", "--output", type=str, help="Prefix output files")
argParser.add_argument("-nf", "--no-filter", action='store_true', help="Filtering by mathcing")
args = argParser.parse_args()

COUNT_NUCLEOTIDE_IN_START = args.count
FILENAMES = args.filenames
OUTPUT = args.output
FILTER_BY_MATCHING=not args.no_filter


def get_complementarity(seq: str) -> str:
    seq = seq.lower()
    replaces = {
        't': 'A',
        'a': 'T',
        'g': 'C',
        'c': 'G'
    }
    for i, j in replaces.items():
        seq = seq.replace(i, j)
    return seq

def get_read_info(filename):
    FILE_SIZE = os.stat(filename).st_size

    _0_5 = {}
    _0_3 = {}
    _16_5 = {}
    _16_3 = {}
    _0_ls = defaultdict(int)
    _16_ls = defaultdict(int)

    progress_step = 0
    all_steps = 302
    step = 0

    with open(filename, 'r') as file:
        for line in file:

            progress_step += len(line)
            if progress_step > FILE_SIZE / all_steps:
                step += progress_step // (FILE_SIZE / all_steps)
                progress_step = progress_step % (FILE_SIZE / all_steps)
                progress(step, all_steps, filename)

            if line[0] == '@':
                continue

            if line.split('	')[2] == '*':
                continue

            if FILTER_BY_MATCHING and not re.match('\d+M$', line.split('	')[5]):
                continue

            seq = line.split('	')[9]
            l = len(seq)

            if line.split('	')[1] == '0':
                _0_ls[l] += 1

                if l in _0_5:
                    _0_5[l] += Counter(seq[:COUNT_NUCLEOTIDE_IN_START])
                else:
                    _0_5[l] = Counter(seq[:COUNT_NUCLEOTIDE_IN_START])

                seq = seq[::-1]
                if l in _0_3:
                    _0_3[l] += Counter(seq[:COUNT_NUCLEOTIDE_IN_START])
                else:
                    _0_3[l] = Counter(seq[:COUNT_NUCLEOTIDE_IN_START])

                continue

            if line.split('	')[1] == '16':
                _16_ls[l] += 1

                seq = get_complementarity(seq[::-1])
                if l in _16_5:
                    _16_5[l] += Counter(seq[:COUNT_NUCLEOTIDE_IN_START])
                else:
                    _16_5[l] = Counter(seq[:COUNT_NUCLEOTIDE_IN_START])

                seq = seq[::-1]
                if l in _16_3:
                    _16_3[l] += Counter(seq[:COUNT_NUCLEOTIDE_IN_START])
                else:
                    _16_3[l] = Counter(seq[:COUNT_NUCLEOTIDE_IN_START])

                continue

        progress(all_steps, all_steps, filename)

    return _0_5, _0_3, _16_5, _16_3, _0_ls, _16_ls

def save_graph(reads_info, graph_name):
    for k in reads_info.keys():
        reads_info[k] = {n: reads_info[k][n]/sum(reads_info[k].values())*100 for n in dict(reads_info[k]).keys()}

    res = [(int(k), v) for k, v in reads_info.items()]
    res = [x for x in res if 18 <= x[0] <= 30]
    res.sort()

    lens_names = [x for x, _ in res]
    nuk_names = ['A', 'C', 'G', 'T']

    fig, ax = plt.subplots()
    bottom = np.zeros(len(lens_names))
    width = 0.9

    for nuk in nuk_names:
        values = [dict(reads_info[l]).get(nuk, 0) for l in lens_names]
        p = ax.bar(
            list(map(str, lens_names)),
            values,
            width,
            label=nuk,
            bottom=bottom
        )
        bottom += values
        ax.bar_label(p, fmt='%.1f', label_type='center', size=7)

    ax.legend(bbox_to_anchor=(1, 1.02), loc="upper left")
    ax.set_title(graph_name)

    plt.subplots_adjust(right=0.87)
    plt.savefig(f'{graph_name}.png')

def save_static_of_lens(reads_info, output_filename):
    # Сохранения количества длин ридов
    res = ""
    for read_info in reads_info:
        filename = read_info[0]
        code_0 = read_info[5]
        code_16 = read_info[6]

        res += filename + '\n'

        res += 'Code 0:\n'
        for i in range(18, 31):
            res += f"{i} - {code_0[i]}\n"
        res += '\n'

        res += 'Code 16:\n'
        for i in range(18, 31):
            res += f"{i} - {code_16[i]}\n"
        res += '\n\n'

    file = open(output_filename, 'w')
    file.write(res)
    file.close()


def to_percents(reads_info):
    reads_info = reads_info.copy()
    for k in reads_info.keys():
        reads_info[k] = {n: reads_info[k][n]/sum(reads_info[k].values())*100 for n in dict(reads_info[k]).keys()}
    return reads_info

if __name__ == '__main__':
    filenames = FILENAMES.split(',')
    if not OUTPUT and len(filenames) > 1:
        raise Exception(
            "При обработке нескольких файлов нужно указать префикс названий выходных файлов. "
            "Укажите через параметр -o. "
            "Пример: -o someprefix"
        )

    if not OUTPUT and len(filenames) == 1:
        OUTPUT = filenames[0]

    reads_info = [(filename, *get_read_info(filename)) for filename in filenames]
    reads_info = [(
        read_info[0],
        to_percents(read_info[1]),
        to_percents(read_info[2]),
        to_percents(read_info[3]),
        to_percents(read_info[4]),
        read_info[5],
        read_info[6],
    ) for read_info in reads_info]
    save_static_of_lens(reads_info, OUTPUT + ".readsinfo.txt")

    _0_5 = {}
    for i in range(18, 31):
        A, T, C, G = 0,0,0,0
        for read_info in reads_info:
            A += read_info[1].get(i, {}).get('A', 0)
            T += read_info[1].get(i, {}).get('T', 0)
            C += read_info[1].get(i, {}).get('C', 0)
            G += read_info[1].get(i, {}).get('G', 0)
        A = A / len(read_info)
        T = T / len(read_info)
        C = C / len(read_info)
        G = G / len(read_info)
        _0_5[i] = {}
        _0_5[i]['A'] = A
        _0_5[i]['T'] = T
        _0_5[i]['C'] = C
        _0_5[i]['G'] = G

    _0_3 = {}
    for i in range(18, 31):
        A, T, C, G = 0,0,0,0
        for read_info in reads_info:
            A += read_info[2].get(i, {}).get('A', 0)
            T += read_info[2].get(i, {}).get('T', 0)
            C += read_info[2].get(i, {}).get('C', 0)
            G += read_info[2].get(i, {}).get('G', 0)
        A = A / len(read_info)
        T = T / len(read_info)
        C = C / len(read_info)
        G = G / len(read_info)
        _0_3[i] = {}
        _0_3[i]['A'] = A
        _0_3[i]['T'] = T
        _0_3[i]['C'] = C
        _0_3[i]['G'] = G

    _16_5 = {}
    for i in range(18, 31):
        A, T, C, G = 0,0,0,0
        for read_info in reads_info:
            A += read_info[3].get(i, {}).get('A', 0)
            T += read_info[3].get(i, {}).get('T', 0)
            C += read_info[3].get(i, {}).get('C', 0)
            G += read_info[3].get(i, {}).get('G', 0)
        A = A / len(read_info)
        T = T / len(read_info)
        C = C / len(read_info)
        G = G / len(read_info)
        _16_5[i] = {}
        _16_5[i]['A'] = A
        _16_5[i]['T'] = T
        _16_5[i]['C'] = C
        _16_5[i]['G'] = G

    _16_3 = {}
    for i in range(18, 31):
        A, T, C, G = 0,0,0,0
        for read_info in reads_info:
            A += read_info[4].get(i, {}).get('A', 0)
            T += read_info[4].get(i, {}).get('T', 0)
            C += read_info[4].get(i, {}).get('C', 0)
            G += read_info[4].get(i, {}).get('G', 0)
        A = A / len(read_info)
        T = T / len(read_info)
        C = C / len(read_info)
        G = G / len(read_info)
        _16_3[i] = {}
        _16_3[i]['A'] = A
        _16_3[i]['T'] = T
        _16_3[i]['C'] = C
        _16_3[i]['G'] = G

    save_graph(_0_5, f'{OUTPUT}.c{COUNT_NUCLEOTIDE_IN_START}.0.5')
    save_graph(_0_3, f'{OUTPUT}.c{COUNT_NUCLEOTIDE_IN_START}.0.3')
    save_graph(_16_5, f'{OUTPUT}.c{COUNT_NUCLEOTIDE_IN_START}.16.5')
    save_graph(_16_3, f'{OUTPUT}.c{COUNT_NUCLEOTIDE_IN_START}.16.3')
