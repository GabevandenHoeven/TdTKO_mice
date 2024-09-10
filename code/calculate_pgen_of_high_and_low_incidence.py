from utils import get_unique_sequences_per_mouse_from_file
from abundance_sequences import get_abundance
import subprocess


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]
    mice_per_file = [13, 10]
    # mice_per_file = [20, 20]

    for file in files:
        max_incidence = mice_per_file[files.index(file)]
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence, seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)

        high_incidence_sequences = seq_per_incidence[-2]
        high_incidence_sequences.extend(seq_per_incidence[-1])
        pgen_input = ['\t'.join([i[1], i[6], i[7]])+'\n' for i in high_incidence_sequences]
        with open('temp.tsv', 'w') as tempfile:
            tempfile.writelines(pgen_input)
        pgen_out = 'pgen-out_high_incidence_' + file.split("\\")[-1]
        command = ['olga-compute_pgen', '-i', 'temp.tsv', '--mouseTRB', '-o', pgen_out, '--v_in', '1', '--j_in', '2']
        subprocess.run(command)

        low_incidence_sequences = seq_per_incidence[0]
        pgen_input = ['\t'.join([i[1], i[6], i[7]]) + '\n' for i in low_incidence_sequences]
        with open('temp.tsv', 'w') as tempfile:
            tempfile.writelines(pgen_input)
        pgen_out = 'pgen-out_low_incidence_' + file.split("\\")[-1]
        command = ['olga-compute_pgen', '-i', 'temp.tsv', '--mouseTRB', '-o', pgen_out, '--v_in', '1', '--j_in', '2']
        subprocess.run(command)

    subprocess.run(['rm', 'temp.tsv'])
