import csv
from plots import plot_d_lengths


def get_d_lengths(filename):
    """
    :param filename:
    :return:
    """
    d_lengths = {
        # D lengths : occurrences
    }
    total = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)
        for line in reader:
            total += 1
            try:
                d_lengths[int(line[header.index('D.length.used')])] += 1
            except KeyError:
                d_lengths.update({int(line[header.index('D.length.used')]): 1})
    return d_lengths, total


if __name__ == '__main__':
    x = []
    y = []
    file = '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv'
    # file = '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv'

    ds, total_lines = get_d_lengths(file)
    x.append(sorted(ds.keys()))
    y.append([ds[e] / total_lines * 100 for e in sorted(ds.keys())])
    file = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
    # file = '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ds, total_lines = get_d_lengths(file)
    x.append(sorted(ds.keys()))
    y.append([ds[e] / total_lines * 100 for e in sorted(ds.keys())])

    out = f'..\\img\\Distribution_of_D_lengths_TdTKO-Normal.png'
    # out = f'..\\img\\Generated_Distribution_of_D_lengths_TdTKO-Normal.png'
    title = f'Distribution of D lengths'
    plot_d_lengths(x, y, ['TdTKO', 'Normal'], title, out)
