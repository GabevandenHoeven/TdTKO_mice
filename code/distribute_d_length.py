import csv
from code.test_scripts.plots import plot_d_lengths


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
    # file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv'
    file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\generated_filtered_data_TdTKO.tsv'
    ds, total_lines = get_d_lengths(file)
    x.append(sorted(ds.keys()))
    y.append([ds[e] / total_lines * 100 for e in sorted(ds.keys())])
    # file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
    file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\generated_filtered_data_Normal.tsv'
    ds, total_lines = get_d_lengths(file)
    x.append(sorted(ds.keys()))
    y.append([ds[e] / total_lines * 100 for e in sorted(ds.keys())])

    # out = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\Distribution_of_D_lengths_TdTKO-Normal.png'
    out = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\' \
          f'Generated_Distribution_of_D_lengths_TdTKO-Normal.png'
    title = f'Distribution of D lengths'
    plot_d_lengths(x, y, ['TdTKO', 'Normal'], title, out)
