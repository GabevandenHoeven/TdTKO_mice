import csv


def get_d_length_and_fraction_with_d(filename):
    """

    :param filename:
    :return:
    """
    avg_d_length = []
    n_lines = 0
    n_d_seg = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)
        for line in reader:
            n_lines += 1
            if int(line[header.index('D.length.used')]) > 0:
                avg_d_length.append(int(line[header.index('D.length.used')]))
                n_d_seg += 1
    fract_d = n_d_seg * 100 / n_lines
    print(f'The average D length is: {sum(avg_d_length) / len(avg_d_length)}')
    print(f'The fraction of sequences that have a D segment is: {n_d_seg} * 100 / {n_lines} = {fract_d}')


if __name__ == '__main__':
    # file = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v1.tsv'
    # file = '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v1.tsv'
    # file = '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v1.tsv'
    # file = '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v1.tsv'
    # file = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
    # file = '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv'
    # file = '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    file = '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv'
    get_d_length_and_fraction_with_d(file)
