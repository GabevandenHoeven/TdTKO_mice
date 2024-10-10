from abundance_sequences import get_abundance
from plots import plot_line_and_scatter_per_incidence
from utils import get_unique_sequences_per_mouse_from_file

if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv'
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO.tsv',
        # '..\\data_files\\WT\\filtered_data\\filtered_data_gen_WT.tsv'
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
                                        '..\\img\\Mean_VJ_distance_per_incidence.png')
    # plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'WT'], ['Fraction of incidence', 'Mean VJ distance (nt)'],
    #                                     'Mean VJ distance per incidence',
    #                                     '..\\img\\Mean_VJ_distance_per_incidence_generated.png')
