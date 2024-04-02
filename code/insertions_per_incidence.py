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

    # for file in files:
    #     x.append([i for i in range(1, mice_per_file[files.index(file)] + 1)])
    #     seqs, n_rows = get_abundance(file, 'd_length >= 0')
    #     insertions = []
    #     for i in range(1, mice_per_file[files.index(file)] + 1):
    #         n_seq_per_incidence = sum(1 for c in seqs.values() if c[0] == i)
    #         try:
    #             mean_ins = sum(c[4][0] for c in seqs.values() if c[0] == i) / n_seq_per_incidence
    #             insertions.append(mean_ins)
    #         except ZeroDivisionError:
    #             mean_ins = 0
    #     y.append(insertions)
    # plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Incidence', 'Mean insertions (nt)'],
    #                                     'Mean insertion length per incidence',
    #                                     'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
    #                                     '\\Mean_insertion_length_per_incidence.png')
    data = []
    for file in files:
        seqs, n_rows = get_abundance(file, 'd_length >= 0')
        insertions = []
        for i in range(1, mice_per_file[files.index(file)] + 1):
            insertions.append([c[4][0] for c in seqs.values() if c[0] == i])
        data.append(insertions)
    plot_boxplot_per_incidence(data[0], ['TdTKO', 'Normal'], ['Incidence', 'Mean insertions (nt)'],
                               'Mean insertion length per incidence',
                               'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                               '\\Boxplot_insertion_length_per_incidence_TdTKO.png')
    plot_boxplot_per_incidence(data[1], ['TdTKO', 'Normal'], ['Incidence', 'Mean insertions (nt)'],
                               'Mean insertion length per incidence',
                               'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                               '\\Boxplot_insertion_length_per_incidence_Normal.png')
