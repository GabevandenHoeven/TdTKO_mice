from utils import get_unique_sequences_from_file
from plots import plot_high_incidence_vj_usage
import numpy


def get_vj_usage_without_d(data):
    """

    :param data:
    :return:
    """
    v_usage = {}
    j_usage = {}
    header = data[0]
    total_lines = 0
    for line in data[1:]:
        v, d_length, j = line[header.index('V.gene')], int(line[header.index('D.length.used')]), \
            line[header.index('J.gene')]
        if d_length == 0:
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
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv',
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]

    v_labels = ['TRBV1*01', 'TRBV12-1*01', 'TRBV12-2*01', 'TRBV13-1*01', 'TRBV13-2*01', 'TRBV13-3*01', 'TRBV14*01',
                'TRBV15*01', 'TRBV16*01', 'TRBV17*01', 'TRBV19*01', 'TRBV2*01', 'TRBV20*01', 'TRBV23*01', 'TRBV24*01',
                'TRBV26*01', 'TRBV29*01', 'TRBV3*01', 'TRBV30*01', 'TRBV31*01', 'TRBV4*01', 'TRBV5*01']
    j_labels = ['TRBJ1-1*01', 'TRBJ1-2*01', 'TRBJ1-3*01', 'TRBJ1-4*01', 'TRBJ1-5*01', 'TRBJ2-1*01', 'TRBJ2-2*01',
                'TRBJ2-3*01', 'TRBJ2-4*01', 'TRBJ2-5*01', 'TRBJ2-7*01']

    v_list = []
    j_list = []
    for file in files:
        filtered_data = get_unique_sequences_from_file(file)
        vs, js, total = get_vj_usage_without_d(filtered_data)
        v_list.append(vs)
        j_list.append(js)
    x_ticks = numpy.arange(1, (len(v_labels)) * 2, 2)
    plot_high_incidence_vj_usage(x_ticks, v_list, (8, 10), 'V Usage in sequences without D segment',
                                 ('V segments', 'Usage (%)'), '..\\img\\V_usage_no_D_segment.png', v_labels)
    x_ticks = numpy.arange(1, (len(j_labels)) * 2, 2)
    plot_high_incidence_vj_usage(x_ticks, j_list, (8, 10), 'J Usage in sequences without D segment',
                                 ('J segments', 'Usage (%)'), '..\\img\\J_usage_no_D_segment.png', j_labels)
