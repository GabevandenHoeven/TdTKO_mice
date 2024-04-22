import csv


def check_sequence_length(filename):
    """

    :param filename:
    :return:
    """
    max_sequence_length = 0
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for line in reader:
            seq_len = len(line[header.index('Junction.nucleotide.sequence')])
            if seq_len > max_sequence_length:
                max_sequence_length = seq_len
    filename = filename.split('\\')[-1]
    print(f'Maximum sequence length for file {filename}:\n{max_sequence_length}')


if __name__ == '__main__':

    file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv'
    check_sequence_length(file)
    file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
    check_sequence_length(file)
