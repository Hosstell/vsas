import json
from collections import Counter
from typing import List, Tuple, Literal

import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

argParser = argparse.ArgumentParser()
argParser.add_argument("-f", "--filename", type=str, required=True, help="Input .sam file")
argParser.add_argument("-c", "--count", type=int, default=1, help="Сount nucleotides in start of a read")
args = argParser.parse_args()

COUNT_NUCLEOTIDE_IN_START = args.count
FILENAME = args.filename

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
    f5 = {}
    f3 = {}
    r5 = {}
    r3 = {}

    with open(filename, 'r') as file:
        for line in file:
            if line[0] == '@':
                continue

            if line.split('	')[2] == '*':
                continue

            seq = line.split('	')[9]
            l = len(seq)

            if line.split('	')[1] == '0':
                if l in f5:
                    f5[l] += Counter(seq[:COUNT_NUCLEOTIDE_IN_START])
                else:
                    f5[l] = Counter(seq[:COUNT_NUCLEOTIDE_IN_START])

                seq = get_complementarity(seq[::-1])
                if l in r3:
                    r3[l] += Counter(seq[:COUNT_NUCLEOTIDE_IN_START])
                else:
                    r3[l] = Counter(seq[:COUNT_NUCLEOTIDE_IN_START])

                continue

            if line.split('	')[1] == '16':
                if l in f3:
                    f3[l] += Counter(seq[:COUNT_NUCLEOTIDE_IN_START])
                else:
                    f3[l] = Counter(seq[:COUNT_NUCLEOTIDE_IN_START])

                seq = get_complementarity(seq[::-1])
                if l in r5:
                    r5[l] += Counter(seq[:COUNT_NUCLEOTIDE_IN_START])
                else:
                    r5[l] = Counter(seq[:COUNT_NUCLEOTIDE_IN_START])

                continue

    return f5, f3, r5, r3

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


if __name__ == '__main__':
    f5, f3, r5, r3 = get_read_info(FILENAME)

    b5 = {i: f5[i] + r5[i] for i in range(18, 31)}
    b3 = {i: f3[i] + r3[i] for i in range(18, 31)}

    save_graph(f5, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.5f')
    save_graph(f3, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.3f')
    save_graph(r5, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.5r')
    save_graph(r3, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.3r')
    save_graph(b5, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.5b')
    save_graph(b3, f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.3b')
