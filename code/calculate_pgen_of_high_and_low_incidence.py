from utils import get_unique_sequences_per_mouse_from_file
from abundance_sequences import get_abundance
import subprocess


def calculate_pgen_for_incidence_group(incidence_sequences, outfile):
    """This function calculates generation probabilities of a given incidence group.
    """
    pgen_input = ['\t'.join([i[1], i[6], i[7]]) + '\n' for i in incidence_sequences]
    with open('temp.tsv', 'w') as tempfile:
        tempfile.writelines(pgen_input)

    command = ['olga-compute_pgen', '-i', 'temp.tsv', '--mouseTRB', '-o', outfile, '--v_in', '1', '--j_in', '2']
    subprocess.run(command)
    return


if __name__ == '__main__':
    files = [
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2 (1).tsv'
    ]

    for file in files:
        max_incidence = 10
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence, seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)

        high_incidence_sequences = seq_per_incidence[-2]
        high_incidence_sequences.extend(seq_per_incidence[-1])
        pgen_out = 'pgen-out_high_incidence_' + file.split("\\")[-1]
        calculate_pgen_for_incidence_group(high_incidence_sequences, pgen_out)

        low_incidence_sequences = seq_per_incidence[0]
        pgen_out = 'pgen-out_low_incidence_' + file.split("\\")[-1]
        calculate_pgen_for_incidence_group(low_incidence_sequences, pgen_out)
