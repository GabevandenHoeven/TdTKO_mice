from utils import get_unique_sequences_per_mouse_from_file
from abundance_sequences import get_abundance
from plots import plot_vj_usage_per_incidence

if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv',
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]
    # segment = 'TRBV1'
    segment = 'TRBV16'
    mice_per_file = [13, 10]
    # There are 13 mice in file 0 from 'files' and 10 in file 1
    # mice_per_file = [20, 20]
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

    title = f'Usage of segment {segment} per incidence for TdTKO and Normal'
    axis_labels = ['Fraction of incidence', 'Percentage of usage (%)']
    out_filename = f'..\\img\\incidence_usage_{segment.replace("*01", "")}'
    plot_vj_usage_per_incidence(x, y, title, axis_labels, out_filename)
