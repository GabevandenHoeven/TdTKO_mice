import csv


def check_sequence_length(filename, delim='\t'):
    """This function checks the sequence length for sequences in a datafile, and prints the length of the longest
    sequence.

    :param filename: str - A filepath to a datafile with CDR3 nucleotide sequences.
    :param delim: str - The delimiter of the datafile. default is tab ('\t')
    :return:
    """
    sequence_lengths = []
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader)
        for line in reader:
            sequence_lengths.append(len(line[header.index('Junction.nucleotide.sequence')]))
        max_sequence_length = max(sequence_lengths)
    filename = filename.split('\\')[-1]
    print(f'Maximum sequence length for file {filename}:\n{max_sequence_length}')


if __name__ == '__main__':

    file = '/data_files/TdTKO/filtered_data/filtered_data_exp_TdTKO_v2.tsv'
    check_sequence_length(file)
    file = '/data_files/WT/filtered_data/filtered_data_exp_WT_v2.tsv'
    check_sequence_length(file)
