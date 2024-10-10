import csv


def count_sequences(filename, delim='\t'):
    """This function counts the number of sequences from each mouse in a datafile.

    :param filename: str - A filepath to a datafile with CDR3 nucleotide sequences.
    :param delim: str - The delimiter of the datafile. default is tab ('\t')
    :return:
    """
    output = {}

    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            mouse = line[header.index('Mouse')]
            try:
                output[mouse] += 1
            except KeyError:
                output.update({mouse: 1})
    return output


if __name__ == '__main__':
    files = [
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv'
    ]
    for file in files:
        out = count_sequences(file)
        print(out)
