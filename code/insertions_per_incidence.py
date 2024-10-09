from abundance_sequences import get_abundance
from plots import plot_line_and_scatter_per_incidence, plot_boxplot_per_incidence
import numpy
from utils import get_unique_sequences_per_mouse_from_file

if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
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
        insertions = []
        for incidence in range(max_incidence):
            incidence_group = seq_per_incidence[incidence]
            try:
                insertions.append(sum(seq[5] for seq in incidence_group) / len(incidence_group))
            except ZeroDivisionError:
                insertions.append(0)
        y.append(insertions)

    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'WT'], ['Fraction of incidence', 'Mean insertions (nt)'],
                                        'Mean insertion length per incidence',
                                        '..\\img\\Mean_insertion_length_per_incidence.png')
    # plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Fraction of incidence', 'Mean insertions (nt)'],
    #                                     'Mean insertion length per incidence',
    #                                     'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
    #                                     '\\Mean_insertion_length_per_incidence_generated.png')

