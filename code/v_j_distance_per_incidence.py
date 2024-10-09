import numpy
from abundance_sequences import get_abundance
from plots import plot_line_and_scatter_per_incidence, plot_boxplot_per_incidence
from utils import get_unique_sequences_per_mouse_from_file

if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal.tsv'
    ]
    mice_per_file = [13, 10]
    # mice_per_file = [20, 20]
    x = []
    y = []

    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence, seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)
        mean_v_j_distances = []
        for i in range(max_incidence):
            incidence_group = seq_per_incidence[i]
            try:
                mean_v_j_distances.append(sum(seq[4] for seq in incidence_group) / len(incidence_group))
            except ZeroDivisionError:
                mean_v_j_distances.append(0)
        y.append(mean_v_j_distances)

    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'WT'], ['Fraction of incidence', 'Mean VJ distance (nt)'],
                                        'Mean VJ distance per incidence',
                                        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                        '\\Mean_VJ_distance_per_incidence.png')
    # plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Fraction of incidence', 'Mean VJ distance (nt)'],
    #                                     'Mean VJ distance per incidence',
    #                                     'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
    #                                     '\\Mean_VJ_distance_per_incidence_generated.png')
    # data = []
    # for data in files:
    #     seqs, n_rows = get_abundance(data, 'd_length >= 0')
    #     vj = []
    #     for i in range(1, mice_per_file[files.index(data)] + 1):
    #         vj.append([c[3][0] for c in seqs.values() if c[0] == i])
    #     data.append(vj)
    # plot_boxplot_per_incidence(data, ['TdTKO', 'Normal'],
    #                            [[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 22, 23], [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]],
    #                            ['Incidence', 'VJ distance (nt)'],
    #                            [1.5, 3.5, 5.5, 7.5, 9.5, 11.5, 13.5, 15.5, 17.5, 19.5, 21, 22, 23],
    #                            list(numpy.arange(1, 14)), 'VJ distance per incidence',
    #                            'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
    #                            '\\Boxplot_VJ_distance_per_incidence.png')

