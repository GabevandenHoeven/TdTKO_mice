import csv


def check_true_d_segment(filename: str, delim='\t'):
    """

    :param filename:
    :param delim:
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
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v4.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v4.tsv',
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v3.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v3.tsv'
    ]
    for f in files:
        fn = f.split('\\')[-1]
        count, total_lines = check_true_d_segment(f)
        fract = count / total_lines
        print(f'The fraction of sequences with correctly identified D segment in file {fn} is: {fract}')
