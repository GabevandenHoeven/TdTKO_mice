from abundance_sequences import get_abundance
from plots import plot_line_and_scatter_per_incidence
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
    v = []
    j = []

    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence, seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)
        v_del = []
        j_del = []
        for incidence in range(max_incidence):
            incidence_group = seq_per_incidence[incidence]
            try:
                v_del.append(sum(seq[8] for seq in incidence_group) / len(incidence_group))
            except ZeroDivisionError:
                v_del.append(0)
            try:
                j_del.append(sum(seq[9] for seq in incidence_group) / len(incidence_group))
            except ZeroDivisionError:
                j_del.append(0)
        v.append(v_del)
        j.append(j_del)

    plot_line_and_scatter_per_incidence(x, v, ['TdTKO', 'Normal'],
                                        ['Fraction of incidence', 'Mean deleted nucleotides (nt)'],
                                        'Mean deleted V gene nucleotides per incidence',
                                        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                        '\\Mean_v_del_per_incidence.png')
    plot_line_and_scatter_per_incidence(x, j, ['TdTKO', 'Normal'],
                                        ['Fraction of incidence', 'Mean deleted nucleotides (nt)'],
                                        'Mean deleted J gene nucleotides per incidence',
                                        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                        '\\Mean_j_del_per_incidence.png')
