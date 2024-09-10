from plots import plot_dist_insertion_length
from utils import get_unique_sequences_from_file
import csv


def get_insertion_counts(file):
    """

    :param fn:
    :param counts:
    :param delim:
    :return:
    """
    total_seq = 0
    insertions = []
    counts = {}
    ins = None
    header = file[0]
    for row in file[1:]:
        total_seq += 1
        try:
            ins = (int(row[header.index('Left.insertion.length')]) - int(row[header.index('Left.palindromic')])) + \
                  (int(row[header.index('Right.insertion.length')]) - int(row[header.index('Right.palindromic')]))
            insertions.append(ins)
            counts[ins] += 1
        except KeyError:
            counts.update({ins: 1})
    return counts, total_seq, insertions


if __name__ == '__main__':
    x_points = []
    y_points = []
    labels = ['TdTKO', 'WT']
    averages = []
    file_list = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv',
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal.tsv'
    ]

    for filename in file_list:
        unique_sequences = get_unique_sequences_from_file(filename)
        insertion_counts, total, insertions = get_insertion_counts(unique_sequences)
        x = sorted(insertion_counts.keys())
        x_points.append(x)
        y_points.append([insertion_counts[point] / total * 100 for point in x])
        averages.append(sum(insertions) / total)
        print(f'The average is: {sum(insertions) / total}')

    outfile = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\' \
            f'Insertion_length_distribution'
    # title = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\' \
    #         f'Generated_Insertion_length_distribution_{"-".join(labels)}'
    plot_dist_insertion_length(x_points, y_points, labels, averages, outfile)
