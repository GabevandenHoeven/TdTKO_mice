from abundance_sequences import get_abundance
from utils import get_unique_sequences_per_mouse_from_file
from plots import number_of_seq_per_incidence

if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv'
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\data_files\\WT\\filtered_data\\filtered_data_gen_WT_v2.tsv'
    ]
    mice_per_file = [13, 10]
    # There are 13 mice in data 0 from 'files' and 10 in data 1
    # mice_per_file = [20, 20]
    x = []
    y = []
    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)[1]
        seq_per_incidence = [len(e) for e in seq_per_incidence]
        total_sequences = len(filtered_data) - 1
        fraction = [seq / total_sequences * 100 for seq in seq_per_incidence]
        x.append([e / max_incidence for e in range(1, max_incidence + 1)])
        y.append(fraction)
    number_of_seq_per_incidence(x, y)
