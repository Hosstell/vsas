import json
from collections import Counter
from typing import List, Tuple, Literal

import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

from progress_bar import progress

argParser = argparse.ArgumentParser()
argParser.add_argument("-f", "--filename", type=str, required=True, help="Input .sam file")
argParser.add_argument("-c", "--count", type=int, default=1, help="Сount nucleotides in start of a read")
args = argParser.parse_args()

COUNT_NUCLEOTIDE_IN_START = args.count
FILENAME = args.filename
FILE_SIZE = os.stat(FILENAME).st_size

def get_read_info(filename):
    _0_5 = {}
    _0_3 = {}
    _16_5 = {}
    _16_3 = {}

    progress_step = 0
    all_steps = 300
    step = 0

    with open(filename, 'r') as file:
        for line in file:

            progress_step += len(line)
            if progress_step > FILE_SIZE / all_steps:
                step += progress_step // (FILE_SIZE / all_steps)
                progress_step = progress_step % (FILE_SIZE / all_steps)
                progress(step, all_steps, FILENAME)

            if line[0] == '@':
                continue

            if line.split('	')[2] == '*':
                continue

            seq = line.split('	')[9]
            l = len(seq)

            if line.split('	')[1] == '0':
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

        progress(all_steps, all_steps, FILENAME)

    return _0_5, _0_3, _16_5, _16_3

def save_graph(reads_info, graph_name):
    for k in reads_info.keys():
        reads_info[k] = {n: reads_info[k][n]/sum(reads_info[k].values())*100 for n in dict(reads_info[k]).keys()}

    res = [(int(k), v) for k, v in reads_info.items()]
    res = [x for x in res if 18 <= x[0] <= 30]
    res.sort()

    file = open(f"{graph_name}.json", 'w')
    file.write(json.dumps(res))
    file.close()

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


if __name__ == '__main__':
    _0_5, _0_3, _16_5, _16_3 = get_read_info(FILENAME)

    save_graph(_0_5, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.0.5')
    save_graph(_0_3, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.0.3')
    save_graph(_16_5, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.16.5')
    save_graph(_16_3, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.16.3')
