from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("-f", "--filename", type=str, required=True, help="Input .sam file")
argParser.add_argument("-c", "--count", type=int, default=1, help="Сount nucleotides in start of a read")
args = argParser.parse_args()

COUNT_NUCLEOTIDE_IN_START = args.count
FILENAME = args.filename


def get_reads_info(filename, count_nucleotide_in_start):
    res_ = {}
    with open(filename, 'r') as file:
        for line in file:
            if line[0] == '@':
                continue

            if line.split('	')[2] == '*':
                continue

            seq = line.split('	')[9]

            len_seq = len(seq)

            if len_seq in res_:
                res_[len_seq] += Counter(seq[:count_nucleotide_in_start])
            else:
                res_[len_seq] = Counter(seq[:count_nucleotide_in_start])
    return res_


def save_graph(read_info, filename):
    for k in read_info.keys():
        read_info[k] = {n: read_info[k][n]/sum(read_info[k].values())*100 for n in dict(read_info[k]).keys()}

    res = [(k, v) for k, v in read_info.items()]
    res = [x for x in res if 13 <= x[0] <= 30]
    res.sort()
    lens_names = [x for x, _ in res]
    nuk_names = ['A', 'T', 'G', 'C', 'N']
    colors = {'A': 'green', 'T': 'red', 'C': 'blue', 'G': 'yellow', 'N': 'black'}

    fig, ax = plt.subplots()
    bottom = np.zeros(len(lens_names))
    width = 0.6

    for nuk in nuk_names:
        values = [dict(read_info[l]).get(nuk, 0) for l in lens_names]
        p = ax.bar(list(map(str, lens_names)), values, width, label=nuk, bottom=bottom)
        bottom += values
        ax.bar_label(p, color=colors[nuk])
        # ax.bar_label(p, label_type='center')

    ax.legend(bbox_to_anchor=(1, 1.02), loc="upper left")

    plt.subplots_adjust(right=0.87)
    plt.savefig(filename)


if __name__ == '__main__':
    reads_info = get_reads_info(FILENAME, COUNT_NUCLEOTIDE_IN_START)
    save_graph(reads_info, f'{FILENAME}.png')
