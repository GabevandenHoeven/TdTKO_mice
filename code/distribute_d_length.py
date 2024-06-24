import csv
from plots import plot_d_lengths
from utils import get_unique_sequences_from_file


def get_d_lengths(file):
    """
    :param file:
    :return:
    """
    d_lengths = {
        # D lengths : occurrences
    }
    total = 0
    header = file[0]

    for line in file[1:]:
        total += 1
        if int(line[header.index('D.length.used')]) not in d_lengths.keys():
            d_lengths.update({int(line[header.index('D.length.used')]): 0})
        d_lengths[int(line[header.index('D.length.used')])] += 1
        # try:
        #     d_lengths[int(line[header.index('D.length.used')])] += 1
        # except KeyError:
        #     d_lengths.update({int(line[header.index('D.length.used')]): 1})
    return d_lengths, total


if __name__ == '__main__':
    x = []
    y = []
    file = '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv'
    # file = '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv'
    filtered_data = get_unique_sequences_from_file(file)
    ds, total_lines = get_d_lengths(filtered_data)
    x.append(sorted(ds.keys()))
    y.append([ds[e] / total_lines * 100 for e in sorted(ds.keys())])
    file = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
    # file = '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    filtered_data = get_unique_sequences_from_file(file)
    ds, total_lines = get_d_lengths(filtered_data)
    x.append(sorted(ds.keys()))
    y.append([ds[e] / total_lines * 100 for e in sorted(ds.keys())])

    out = f'..\\img\\unique_seq_img\\Distribution_of_D_lengths_TdTKO-Normal.png'
    # out = f'..\\img\\unique_seq_img\\Generated_Distribution_of_D_lengths_TdTKO-Normal.png'
    title = f'Distribution of D lengths'
    plot_d_lengths(x, y, ['TdTKO', 'Normal'], title, out)
