import json
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

argParser = argparse.ArgumentParser()
argParser.add_argument("-f", "--filename", type=str, required=True, help="Input .sam file")
argParser.add_argument("-c", "--count", type=int, default=1, help="Сount nucleotides in start of a read")
argParser.add_argument("-r", "--reverse", action='store_true', help="Show reverse?")
args = argParser.parse_args()

COUNT_NUCLEOTIDE_IN_START = args.count
FILENAME = args.filename
IS_REVERSE = args.reverse

def get_reads_info_forward(filename, count_nucleotide_in_start):
    res_ = {}
    with open(filename, 'r') as file:
        for line in file:
            if line[0] == '@':
                continue

            if line.split('	')[2] == '*':
                continue

            if line.split('	')[1] == '16':
                continue

            seq = line.split('	')[9]
            len_seq = len(seq)

            if len_seq in res_:
                res_[len_seq] += Counter(seq[:count_nucleotide_in_start])
            else:
                res_[len_seq] = Counter(seq[:count_nucleotide_in_start])
    return res_


def get_reads_info_reverse(filename, count_nucleotide_in_start):
    res_ = {}
    with open(filename, 'r') as file:
        for line in file:
            if line[0] == '@':
                continue

            if line.split('	')[2] == '*':
                continue

            if line.split('	')[1] == '0':
                continue

            seq = line.split('	')[9]

            seq = seq[::-1].lower()
            replaces = {
                't': 'A',
                'a': 'T',
                'g': 'C',
                'c': 'G'
            }
            for i, j in replaces.items():
                seq = seq.replace(i, j)

            len_seq = len(seq)

            if len_seq in res_:
                res_[len_seq] += Counter(seq[:count_nucleotide_in_start])
            else:
                res_[len_seq] = Counter(seq[:count_nucleotide_in_start])
    return res_

def get_reads_info(filename, count_nucleotide_in_start, is_reverse):
    if is_reverse:
        return get_reads_info_reverse(filename, count_nucleotide_in_start)
    else:
        return get_reads_info_forward(filename, count_nucleotide_in_start)

def get_reads_info_debug():
    devdata = json.loads(open("devdata.json", "r").read())
    return {int(k): v for k, v in devdata.items()}


def save_graph(read_info, filename, direction):
    for k in read_info.keys():
        read_info[k] = {n: read_info[k][n]/sum(read_info[k].values())*100 for n in dict(read_info[k]).keys()}

    res = [(int(k), v) for k, v in read_info.items()]
    res = [x for x in res if 13 <= x[0] <= 30]
    res.sort()
    lens_names = [x for x, _ in res]
    nuk_names = ['A', 'C', 'G', 'T']

    fig, ax = plt.subplots()
    bottom = np.zeros(len(lens_names))
    width = 0.9

    for nuk in nuk_names:
        values = [dict(read_info[l]).get(nuk, 0) for l in lens_names]
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
    ax.set_title(f'{FILENAME} (-c {COUNT_NUCLEOTIDE_IN_START}, {direction})')

    plt.subplots_adjust(right=0.87)
    plt.savefig(filename)


if __name__ == '__main__':
    direction = 'reverse' if IS_REVERSE else 'forward'
    filename = f'{FILENAME}.c{COUNT_NUCLEOTIDE_IN_START}.{direction}.png'

    if os.getenv('DEV'):
        reads_info = get_reads_info_debug()
    else:
        reads_info = get_reads_info(FILENAME, COUNT_NUCLEOTIDE_IN_START, IS_REVERSE)
    save_graph(reads_info, filename, direction)
