from utils import get_unique_sequences_per_mouse_from_file
from abundance_sequences import get_abundance
from plots import plot_vj_usage_per_incidence


def get_incidence_segment_usage(segment, files, mice_per_file):
    """Plots the incidence usage of a given V or J gene segment.

    :param segment: str - The segment for which we want to plot the usage over incidence.
    :param files: list - A list of filepaths to the datafiles.
    :param mice_per_file: list - A list of the number of mice in each of the corresponding datafiles from `files`.
    :return:
    """
    x = []
    y = []
    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        x.append([i for i in range(1, max_incidence + 1)])
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)[1]
        usage_values = []
        for incidence_group in seq_per_incidence:
            v_segment_count = sum([1 for seq in incidence_group if seq[6] == segment])
            usage_values.append(v_segment_count / len(incidence_group) * 100)
        y.append(usage_values)

    segment = segment.replace('*01', '')
    title = f'Usage of segment {segment} per incidence for TdTKO and WT'
    axis_labels = ['Fraction of incidence', 'Percentage of usage (%)']
    out_filename = f'..\\img\\incidence_usage_{segment.replace("*01", "")}'
    plot_vj_usage_per_incidence(x, y, title, axis_labels, out_filename)
    return


if __name__ == '__main__':
    fs = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv',
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\data_files\\WT\\filtered_data\\filtered_data_gen_WT_v2.tsv'
    ]
    mpf = [13, 10]
    # There are 13 mice in data 0 from 'files' and 10 in data 1

    s = 'TRBV1*01'
    get_incidence_segment_usage(s, fs, mpf)
    s = 'TRBV16*01'
    get_incidence_segment_usage(s, fs, mpf)

