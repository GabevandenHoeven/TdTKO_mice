import csv


def check_true_d_segment(filename: str, delim='\t'):
    """This function opens a datafile of synthetic sequences and compares the inferred D segment length to the true
    D segment length

    :param filename: str - The filepath of the datafile.
    :param delim: str - The delimiter of the datafile. Default is tab ('\t')
    :return:
    """
    count_correct_d_segment = 0
    lines = 0
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=delim)
        header = next(reader)
        for row in reader:
            lines += 1
            called_d, true_d = row[header.index('D.used')], row[header.index('D.used.data')]
            if called_d == true_d:
                count_correct_d_segment += 1
    return count_correct_d_segment, lines


if __name__ == '__main__':
    files = [
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v1.tsv',
        '..\\..\\data_files\\WT\\filtered_data\\filtered_data_gen_WT_v1.tsv',
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        '..\\..\\data_files\\WT\\filtered_data\\filtered_data_gen_WT_v2.tsv'
    ]
    for f in files:
        fn = f.split('\\')[-1]
        count, total_lines = check_true_d_segment(f)
        fract = count / total_lines
        print(f'The fraction of sequences with correctly identified D segment in file {fn} is: {fract}')
