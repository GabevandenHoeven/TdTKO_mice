from abundance_sequences import get_abundance
from utils import get_unique_sequences_per_mouse_from_file
from plots import plot_high_incidence_vj_usage
import numpy

if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal.tsv'
    ]
    mice_per_file = [13, 10]
    # mice_per_file = [20, 20]
    v_list = []
    j_list = []

    v_labels = ['TRBV1*01', 'TRBV12-1*01', 'TRBV12-2*01', 'TRBV13-1*01', 'TRBV13-2*01', 'TRBV13-3*01', 'TRBV14*01',
                'TRBV15*01', 'TRBV16*01', 'TRBV17*01', 'TRBV19*01', 'TRBV2*01', 'TRBV20*01', 'TRBV23*01', 'TRBV24*01',
                'TRBV26*01', 'TRBV29*01', 'TRBV3*01', 'TRBV30*01', 'TRBV31*01', 'TRBV4*01', 'TRBV5*01']
    j_labels = ['TRBJ1-1*01', 'TRBJ1-2*01', 'TRBJ1-3*01', 'TRBJ1-4*01', 'TRBJ1-5*01', 'TRBJ2-1*01', 'TRBJ2-2*01',
                'TRBJ2-3*01', 'TRBJ2-4*01', 'TRBJ2-5*01', 'TRBJ2-7*01']

    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence, seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)
        high_incidence_sequences = seq_per_incidence[-2]
        high_incidence_sequences.extend(seq_per_incidence[-1])
        # high_incidence_sequences = seq_per_incidence[0] # low incidence, one mice
        vs = [seq[6] for seq in high_incidence_sequences]
        js = [seq[7] for seq in high_incidence_sequences]
        v = {}
        j = {}
        for i in set(vs):
            v.update({i: sum(1 for ii in vs if ii == i) / len(high_incidence_sequences) * 100})
        for i in set(js):
            j.update({i: sum(1 for ii in js if ii == i) / len(high_incidence_sequences) * 100})
        for i in v_labels:
            if i not in v.keys():
                v.update({i: 0})
        for i in j_labels:
            if i not in j.keys():
                j.update({i: 0})
        v_list.append([v[e] for e in sorted(v.keys())])
        j_list.append([j[e] for e in sorted(j.keys())])

    x_ticks = numpy.arange(1, (len(v_labels)) * 2, 2)
    plot_high_incidence_vj_usage(x_ticks, v_list, (8, 10), 'V Usage High Incidence',
                                 ('V segments', 'Usage (%)'), '..\\img\\V_usage_high_incidence.png', v_labels)
    # plot_high_incidence_vj_usage(x_ticks, v_list, (8, 10), 'V Usage Low Incidence',
    #                              ('V segments', 'Usage (%)'), '..\\img\\V_usage_low_incidence.png', v_labels)
    #
    x_ticks = numpy.arange(1, (len(j_labels)) * 2, 2)
    plot_high_incidence_vj_usage(x_ticks, j_list, (8, 10), 'J Usage High Incidence',
                                 ('J segments', 'Usage (%)'), '..\\img\\J_usage_high_incidence.png', j_labels)
    # plot_high_incidence_vj_usage(x_ticks, j_list, (8, 10), 'J Usage Low Incidence',
    #                              ('J segments', 'Usage (%)'), '..\\img\\J_usage_low_incidence.png', j_labels)
