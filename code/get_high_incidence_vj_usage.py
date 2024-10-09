from abundance_sequences import get_abundance
from utils import get_unique_sequences_per_mouse_from_file
from plots import plot_high_incidence_vj_usage
import numpy


def incidence_vj_usage(incidence_sequences, segment_labels, index, plot_title, axes_labels, output_filepath):
    """This function takes a group of sequences from a certain incidence and plots the V or J segment usage of these
    sequences as well as showing the fraction of sequences that have no inferred D segment.

    :param incidence_sequences: list - A list of sequence information for sequences of a certain incidence.
    :param segment_labels: list - A list of segments to be shown in the plot
    :param index: int - the corresponding index to the segments
    :param plot_title: str - The title of the plot.
    :param axes_labels: list/tuple - the labels for the x-axis (0) and y-axis (1).
    :param output_filepath: str - The filepath to the output image.

    :return:
    """
    segments_list = []
    no_d_for_segment = []
    segments = [seq[index] for seq in incidence_sequences]

    d_lengths_segments = {}
    for seq in incidence_sequences:
        try:
            d_lengths_segments[seq[index].replace('*01', '')][seq[3]] += 1
        except KeyError:
            d_lengths_segments.update({seq[index].replace('*01', ''): [0] * 15})
            d_lengths_segments[seq[index].replace('*01', '')][seq[3]] += 1

    seg = {}
    for i in set(segments):
        seg.update({i.replace('*01', ''): sum(1 for ii in segments if ii == i) / len(incidence_sequences) * 100})
    for i in segment_labels:
        if i not in seg.keys():
            seg.update({i: 0})
        if i not in d_lengths_segments.keys():
            d_lengths_segments.update({i: [0] * 15})

    segments_list.append([seg[e] for e in sorted(seg.keys())])
    no_d_for_segment.append([d_lengths_segments[i][0] / len(high_incidence_sequences) * 100
                             for i in sorted(d_lengths_segments.keys())])

    x_ticks = numpy.arange(1, (len(segment_labels)) * 2, 2)
    plot_high_incidence_vj_usage(x_ticks, segments_list, no_d_for_segment, (8, 10),
                                 plot_title, axes_labels, output_filepath, segment_labels)
    return


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]
    mice_per_file = [13, 10]
    # mice_per_file = [20, 20]

    v_segments = ['TRBV1*01', 'TRBV12-1*01', 'TRBV12-2*01', 'TRBV13-1*01', 'TRBV13-2*01', 'TRBV13-3*01', 'TRBV14*01',
                'TRBV15*01', 'TRBV16*01', 'TRBV17*01', 'TRBV19*01', 'TRBV2*01', 'TRBV20*01', 'TRBV23*01', 'TRBV24*01',
                'TRBV26*01', 'TRBV29*01', 'TRBV3*01', 'TRBV30*01', 'TRBV31*01', 'TRBV4*01', 'TRBV5*01']

    j_segments = ['TRBJ1-1*01', 'TRBJ1-2*01', 'TRBJ1-3*01', 'TRBJ1-4*01', 'TRBJ1-5*01', 'TRBJ2-1*01', 'TRBJ2-2*01',
                'TRBJ2-3*01', 'TRBJ2-4*01', 'TRBJ2-5*01', 'TRBJ2-7*01']

    # The following are all the segments including the non-functional segments.
    # v_segments = [
    #     'TRBV1', 'TRBV10', 'TRBV11', 'TRBV12-1', 'TRBV12-2', 'TRBV12-3', 'TRBV13-1', 'TRBV13-2', 'TRBV13-3', 'TRBV14',
    #     'TRBV15', 'TRBV16', 'TRBV17', 'TRBV18', 'TRBV19', 'TRBV2', 'TRBV20', 'TRBV21', 'TRBV22', 'TRBV23', 'TRBV24',
    #     'TRBV25', 'TRBV26', 'TRBV27', 'TRBV28', 'TRBV29', 'TRBV3', 'TRBV30', 'TRBV31', 'TRBV4', 'TRBV5', 'TRBV6',
    #     'TRBV7', 'TRBV8', 'TRBV9']
    # j_segments = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ1-6', 'TRBJ1-7',
    #             'TRBJ2-1', 'TRBJ2-2', 'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-6', 'TRBJ2-7']

    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence, seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)
        high_incidence_sequences = seq_per_incidence[-2]
        high_incidence_sequences.extend(seq_per_incidence[-1])
        low_incidence_sequences = seq_per_incidence[0]

        title = 'V Usage High Incidence'
        labels = ('V segments', 'Usage (%)')
        outfile = '..\\img\\V_usage_high_incidence.png'
        incidence_vj_usage(high_incidence_sequences, v_segments, 6, title, labels, outfile)
        title = 'J Usage High Incidence'
        labels = ('J segments', 'Usage (%)')
        outfile = '..\\img\\J_usage_high_incidence.png'
        incidence_vj_usage(high_incidence_sequences, j_segments, 7, title, labels, outfile)

        title = 'V Usage Low Incidence'
        labels = ('V segments', 'Usage (%)')
        outfile = '..\\img\\V_usage_low_incidence.png'
        incidence_vj_usage(low_incidence_sequences, v_segments, 6, title, labels, outfile)
        title = 'J Usage Low Incidence'
        labels = ('J segments', 'Usage (%)')
        outfile = '..\\img\\J_usage_low_incidence.png'
        incidence_vj_usage(low_incidence_sequences, j_segments, 7, title, labels, outfile)

