from utils import get_unique_sequences_from_file
from plots import plot_vj_usage_with_without_d
import numpy


def get_vj_usage_with_or_without_d(data, exp: str):
    """This function get the VJ usage for sequences with an inferred D length based on the given expression.

    :param data: nested list - A nested list of lines from a datafile. The first line is a header.
    :param exp: A string with a boolean statement to evaluate.
    :return:
    """
    v_usage = {}
    j_usage = {}
    header = data[0]
    total_lines = 0
    for line in data[1:]:
        v, d_length, j = line[header.index('V.gene')], int(line[header.index('D.length.used')]), \
            line[header.index('J.gene')]
        if eval(exp):
            total_lines += 1
            try:
                v_usage[v] += 1
            except KeyError:
                v_usage.update({v: 1})
            try:
                j_usage[j] += 1
            except KeyError:
                j_usage.update({j: 1})
    v_usage = [v_usage[v] / total_lines * 100 for v in sorted(v_usage.keys())]
    j_usage = [j_usage[j] / total_lines * 100 for j in sorted(j_usage.keys())]
    return v_usage, j_usage, total_lines


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv',
    ]

    v_labels = ['TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV13-1', 'TRBV13-2', 'TRBV13-3', 'TRBV14',
                'TRBV15', 'TRBV16', 'TRBV17', 'TRBV19', 'TRBV2', 'TRBV20', 'TRBV23', 'TRBV24',
                'TRBV26', 'TRBV29', 'TRBV3', 'TRBV30', 'TRBV31', 'TRBV4', 'TRBV5']
    j_labels = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ2-1', 'TRBJ2-2',
                'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-7']

    v_list = []
    j_list = []
    for file in files:
        filtered_data = get_unique_sequences_from_file(file)
        vs, js, total = get_vj_usage_with_or_without_d(filtered_data, 'd_length == 0')
        v_list.append(vs)
        j_list.append(js)
        print(f'Total number of unique sequences in all mice combined without D: {total}')
        vs, js, total = get_vj_usage_with_or_without_d(filtered_data, 'd_length != 0')
        v_list.append(vs)
        j_list.append(js)
        print(f'Total number of unique sequences in all mice combined with D: {total}')
    x_ticks = numpy.arange(1, (len(v_labels)) * 2, 2)
    plot_vj_usage_with_without_d(x_ticks, [v_list[2], v_list[3]], (8, 10),
                                 'V Usage in WT sequences with and without D segment',
                                 ('V segments', 'Usage (%)'), '..\\img\\V_usage_with_without_D_segment.png', v_labels)
    x_ticks = numpy.arange(1, (len(j_labels)) * 2, 2)
    plot_vj_usage_with_without_d(x_ticks, [j_list[2], j_list[3]], (8, 10),
                                 'J Usage in WT sequences with and without D segment',
                                 ('J segments', 'Usage (%)'), '..\\img\\J_usage_with_without_D_segment.png', j_labels)
