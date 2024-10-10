from abundance_sequences import get_abundance
from plots import plot_line_and_scatter_per_incidence
from utils import get_unique_sequences_per_mouse_from_file


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv'
    ]
    mice_per_file = [13, 10]
    # There are 13 mice in data 0 from 'files' and 10 in data 1

    x = []
    y = []

    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)[1]
        mean_d_lengths = []
        for incidence in range(max_incidence):
            incidence_group = seq_per_incidence[incidence]
            try:
                mean_d_lengths.append(sum(seq[3] for seq in incidence_group) / len(incidence_group))
            except ZeroDivisionError:
                mean_d_lengths.append(0)
        y.append(mean_d_lengths)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'WT'],
                                        ['Fraction of incidence', 'Inferred D-segment length (nt)'],
                                        'Mean inferred D-length per incidence',
                                        '..\\img\\Mean_inferred_D_length_per_incidence.png')

    x = []
    y = []
    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence = get_abundance(filtered_data, 'seq[3] == 0', max_incidence)[0]
        y.append(fract_incidence)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'WT'],
                                        ['Fraction of incidence', 'Percentage of sequences (%)'],
                                        'Inferred D-length = 0 nt',
                                        '..\\img\\Fraction_no_D_per_incidence.png')

    x = []
    y = []
    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence = get_abundance(filtered_data, 'seq[3] <= 2', max_incidence)[0]
        y.append(fract_incidence)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'WT'],
                                        ['Fraction of incidence', 'Percentage of sequences (%)'],
                                        'Inferred D-length <= 2 nt',
                                        '..\\img\\Fraction_max_2nt_D_per_incidence.png')
