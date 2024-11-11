import csv
import subprocess


def compute_pgen_for_data_file(filename, out_filename, delim='\t'):
    """

    :param filename:
    :param out_filename:
    :param delim:
    :return:
    """
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        sequence_index = header.index('nucleotide_sequence')
    command = ['olga-compute_pgen', '-i', filename, '--lines_to_skip', '1', '--humanTRB', '-o', out_filename,
               '--seq_in', str(sequence_index)]
    subprocess.run(command)
    return


if __name__ == '__main__':
    fn1 = '..\\data_files\\Tregs_data\\control\\CD4_Treg\\filtered_data\\filtered_data_Tregs.tsv'
    out_fn1 = '..\\data_files\\Tregs_data\\control\\CD4_Treg\\filtered_data\\pgen_output.tsv'
    compute_pgen_for_data_file(fn1, out_fn1)

    fn2 = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\filtered_data_TNaive.tsv'
    out_fn2 = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\pgen_output.tsv'
    # out_fn2 = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\pgen_output (1).tsv'
    # out_fn2 = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\pgen_output (2).tsv'
    compute_pgen_for_data_file(fn2, out_fn2)
