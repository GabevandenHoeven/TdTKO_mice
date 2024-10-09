import csv


def check_insertions(filename: str, delim='\t'):
    """Checks if the sequences in the data contain insertions.
    :param filename:
    :param delim: str
    :return:
    """
    insertion_count = 0
    line_count = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            line_count += 1
            l_ins, r_ins, l_palindromic, r_palindromic = int(line[header.index('Left.insertion.length')]), \
                int(line[header.index('Right.insertion.length')]), int(line[header.index('Left.palindromic')]), \
                int(line[header.index('Right.palindromic')])
            if (l_ins - l_palindromic) + (r_ins - r_palindromic) != 0:
                insertion_count += 1

    return insertion_count, line_count


if __name__ == '__main__':
    files = [
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v1.tsv',
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v1.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv',
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v1.tsv',
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v1.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]
    for file in files:
        insertions, lines = check_insertions(file)
        fn = file.split('\\')[-1]
        print(f'Number of sequences with insertions in data \"{fn}\": {insertions}. {insertions / lines}')
