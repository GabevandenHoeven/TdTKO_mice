import numpy
from abundance_sequences import get_abundance
from plots import plot_line_and_scatter_per_incidence, plot_boxplot_per_incidence


if __name__ == '__main__':
    files = [
        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv',
        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
    ]
    mice_per_file = [13, 10]
    x = []
    y = []

    for file in files:
        x.append([i for i in range(1, mice_per_file[files.index(file)] + 1)])
        seqs, n_rows = get_abundance(file, 'd_length >= 0')
        mean_v_j_distances = []
        n_seq = []
        for i in range(1, mice_per_file[files.index(file)] + 1):
            n_seq_per_incidence = sum(1 for c in seqs.values() if c[0] == i)
            n_seq.append(n_seq_per_incidence)
            try:
                mean_vj = sum(c[3][0] for c in seqs.values() if c[0] == i) / n_seq_per_incidence
                mean_v_j_distances.append(mean_vj)
            except ZeroDivisionError:
                mean_vj = 0
        print(n_seq)
        y.append(mean_v_j_distances)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Incidence', 'Mean VJ distance (nt)'],
                                        'Mean VJ distance per incidence',
                                        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                        '\\Mean_VJ_distance_per_incidence.png')
    data = []
    for file in files:
        seqs, n_rows = get_abundance(file, 'd_length >= 0')
        vj = []
        for i in range(1, mice_per_file[files.index(file)] + 1):
            vj.append([c[3][0] for c in seqs.values() if c[0] == i])
        data.append(vj)
    plot_boxplot_per_incidence(data, ['TdTKO', 'Normal'],
                               [[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 22, 23], [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]],
                               ['Incidence', 'VJ distance (nt)'],
                               [1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21, 22, 23],
                               list(numpy.arange(1, 14)), 'VJ distance per incidence',
                               'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                               '\\Boxplot_VJ_distance_per_incidence.png')

