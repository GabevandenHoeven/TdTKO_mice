from abundance_sequences import get_abundance
from plots import plot_line_and_scatter_per_incidence, plot_boxplot_per_incidence
import numpy

if __name__ == '__main__':
    files = [
        # 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv',
        # 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
        '..\\data_files\\TdTKo\\filtered_data\\generated_filtered_data_TdTKO.tsv',
        '..\\data_files\\B6\\filtered_data\\generated_filtered_data_Normal.tsv'
    ]
    # mice_per_file = [13, 10]
    # There are 13 mice in file 0 from 'files' and 10 in file 1
    mice_per_file = [20, 20]

    x = []
    y = []

    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        seq_per_incidence = get_abundance(file, 'seq[3] >= 0', max_incidence)[1]
        mean_d_lengths = []
        for incidence in range(max_incidence):
            incidence_group = seq_per_incidence[incidence]
            try:
                mean_d_lengths.append(sum(seq[3] for seq in incidence_group) / len(incidence_group))
            except ZeroDivisionError:
                mean_d_lengths.append(0)
        y.append(mean_d_lengths)

    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Incidence', 'Inferred D-segment length (nt)'],
                                        'Mean inferred D-length',
                                        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                        '\\Mean_inferred_D_length_per_incidence_generated.png')

    x = []
    y = []
    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        fract_incidence = get_abundance(file, 'seq[3] == 0', max_incidence)[0]
        y.append(fract_incidence)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Incidence', 'Percentage of sequences (%)'],
                                        'Inferred D-length = 0 nt',
                                        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                        '\\Fraction_no_D_per_incidence_generated.png')
    x = []
    y = []
    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        fract_incidence = get_abundance(file, 'seq[3] <= 2', max_incidence)[0]
        y.append(fract_incidence)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Incidence', 'Percentage of sequences (%)'],
                                        'Inferred D-length <= 2 nt',
                                        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                        '\\Fraction_max_2nt_D_per_incidence_generated.png')

    # data = []
    # for file in files:
    #     seqs, n_rows = get_abundance(file, 'd_length >= 0')
    #     d_lengths = []
    #     for i in range(1, mice_per_file[files.index(file)] + 1):
    #         d_lengths.append([c[1][0] for c in seqs.values() if c[0] == i])
    #     data.append(d_lengths)
    # plot_boxplot_per_incidence(data, ['TdTKO', 'Normal'],
    #                            [[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 22, 23], [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]],
    #                            ['Incidence', 'D length (nt)'],
    #                            [1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21, 22, 23],
    #                            list(numpy.arange(1, 14)), 'D length per incidence',
    #                            'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
    #                            '\\Boxplot_D_length_per_incidence.png')
