from plots import plot_fractions_ins
from utils import get_unique_sequences_from_file


def get_fraction_per_read_threshold(file: str, thresholds: list):
    """Calculates the fraction of sequences found containing insertions having above a given threshold of supporting
    reads. Returns a list of fractions, one for each of the thresholds in the parameter.

    :param file: str - A data path to a datafile.
    :param thresholds: list - A list containing thresholds for a minimum number of reads.
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


def get_fraction_per_read_threshold_old(file: str, thresholds: list):
    """Calculates the fraction of sequences found containing insertions having above a given threshold of supporting
    reads. Returns a list of fractions, one for each of the thresholds in the parameter. This function does not account
    for P nucleotides.

    :param file: str - A data path to a datafile.
    :param thresholds: list - A list containing thresholds for a minimum number of reads.
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
                ins = int(line[header.index('Insertions')])
                if ins > 0:
                    count_ins += 1
        fraction = count_ins / count_lines_supporting_reads * 100
        y.append(fraction)
    return y


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv',
    ]
    thres = [1, 2, 3, 5, 10, 15, 20, 25, 50]
    y_points = []
    labels = ['TdTKO', 'WT']
    for file in files:
        unique_sequences = get_unique_sequences_from_file(file)
        fracts = get_fraction_per_read_threshold(unique_sequences, thres)
        y_points.append(fracts)

    title = '..\\img\\Fractions_insertions_per_reads.png'

    plot_fractions_ins(thres, y_points, title, labels)

    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v1.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v1.tsv',
    ]
    y_points = []
    for file in files:
        unique_sequences = get_unique_sequences_from_file(file)
        fracts = get_fraction_per_read_threshold_old(unique_sequences, thres)
        y_points.append(fracts)
    title = '..\\img\\Fractions_insertions_per_reads_old.png'
    plot_fractions_ins(thres, y_points, title, labels)
