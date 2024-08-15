from plots import plot_fractions_ins
import csv
from utils import get_unique_sequences_from_file


def get_fraction_per_read_threshold(file: str, thresholds: list, delim='\t'):
    """

    :param filename:
    :param thresholds:
    :param delim:
    :return:
    """
    y = []
    for threshold in thresholds:
        count_lines_supporting_reads = 0
        count_ins = 0
        header = file[0]
        for line in file[1:]:
            if int(line[header.index('Number.of.reads')]) > threshold:
                count_lines_supporting_reads += 1
                ins = (int(line[header.index('Left.insertion.length')]) -
                       int(line[header.index('Left.palindromic')])) + \
                      (int(line[header.index('Right.insertion.length')]) -
                       int(line[header.index('Right.palindromic')]))
                if ins > 0:
                    count_ins += 1
        fraction = count_ins / count_lines_supporting_reads * 100
        y.append(fraction)
    return y


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv',
    ]
    thres = [0, 1, 2, 5, 10, 15, 20, 25, 50]
    y_points = []
    for file in files:
        unique_sequences = get_unique_sequences_from_file(file)
        fracts = get_fraction_per_read_threshold(unique_sequences, thres)
        y_points.append(fracts)

    title = '..\\img\\Normal-TdTKO_fraction_insertions_per_reads (1).png'
    plot_fractions_ins(thres, y_points, title)
