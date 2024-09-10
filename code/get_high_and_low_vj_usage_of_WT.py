from abundance_sequences import get_abundance
from utils import get_unique_sequences_per_mouse_from_file
from plots import plot_high_low_vj_usage
import numpy

if __name__ == '__main__':
    file = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'

    mice_in_file = 10
    v_list = []
    j_list = []
    no_d_for_v = []
    no_d_for_j = []

    v_labels = ['TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV13-1', 'TRBV13-2', 'TRBV13-3', 'TRBV14',
                'TRBV15', 'TRBV16', 'TRBV17', 'TRBV19', 'TRBV2', 'TRBV20', 'TRBV23', 'TRBV24',
                'TRBV26', 'TRBV29', 'TRBV3', 'TRBV30', 'TRBV31', 'TRBV4', 'TRBV5']
    # v_labels = [
    #     'TRBV1', 'TRBV10', 'TRBV11', 'TRBV12-1', 'TRBV12-2', 'TRBV12-3', 'TRBV13-1', 'TRBV13-2', 'TRBV13-3', 'TRBV14',
    #     'TRBV15', 'TRBV16', 'TRBV17', 'TRBV18', 'TRBV19', 'TRBV2', 'TRBV20', 'TRBV21', 'TRBV22', 'TRBV23', 'TRBV24',
    #     'TRBV25', 'TRBV26', 'TRBV27', 'TRBV28', 'TRBV29', 'TRBV3', 'TRBV30', 'TRBV31', 'TRBV4', 'TRBV5', 'TRBV6',
    #     'TRBV7', 'TRBV8', 'TRBV9']
    j_labels = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ2-1', 'TRBJ2-2',
                'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-7']
    # j_labels = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ1-7',
    #             'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7']

    max_incidence = mice_in_file
    filtered_data = get_unique_sequences_per_mouse_from_file(file)
    fract_incidence, seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)

    high_incidence_sequences = seq_per_incidence[-2]
    high_incidence_sequences.extend(seq_per_incidence[-1])
    vs = [seq[6] for seq in high_incidence_sequences]
    js = [seq[7] for seq in high_incidence_sequences]
    v = {}
    j = {}
    for i in set(vs):
        v.update({i.replace('*01', ''): sum(1 for ii in vs if ii == i) / len(high_incidence_sequences) * 100})
    for i in v_labels:
        if i not in v.keys():
            v.update({i: 0})
    for i in set(js):
        j.update({i.replace('*01', ''): sum(1 for ii in js if ii == i) / len(high_incidence_sequences) * 100})
    for i in j_labels:
        if i not in j.keys():
            j.update({i: 0})
    v_list.append([v[e] for e in sorted(v.keys())])
    j_list.append([j[e] for e in sorted(j.keys())])

    low_incidence_sequences = seq_per_incidence[0]  # low incidence, one mouse
    vs = [seq[6] for seq in low_incidence_sequences]
    js = [seq[7] for seq in low_incidence_sequences]
    v = {}
    j = {}
    for i in set(vs):
        v.update({i.replace('*01', ''): sum(1 for ii in vs if ii == i) / len(low_incidence_sequences) * 100})
    for i in v_labels:
        if i not in v.keys():
            v.update({i: 0})
    for i in set(js):
        j.update({i.replace('*01', ''): sum(1 for ii in js if ii == i) / len(low_incidence_sequences) * 100})
    for i in j_labels:
        if i not in j.keys():
            j.update({i: 0})
    v_list.append([v[e] for e in sorted(v.keys())])
    j_list.append([j[e] for e in sorted(j.keys())])

    x_ticks = numpy.arange(1, (len(v_labels)) * 2, 2)
    plot_high_low_vj_usage(x_ticks, v_list, (8, 10), 'V Usage High Incidence vs Low Incidence\nIn WT',
                           ('V segments', 'Usage (%)'), '..\\img\\V_usage_high_and_low_incidence.png', v_labels)

    x_ticks = numpy.arange(1, (len(j_labels)) * 2, 2)
    plot_high_low_vj_usage(x_ticks, j_list, (8, 10), 'J Usage High Incidence vs Low Incidence\nIn WT',
                           ('J segments', 'Usage (%)'), '..\\img\\J_usage_high_and_low_incidence.png', j_labels)
