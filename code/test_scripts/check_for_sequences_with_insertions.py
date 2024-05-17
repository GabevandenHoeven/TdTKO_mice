import csv


def check_insertions(filename: str, delim='\t'):
    """Checks if the sequences in the file contain insertions.
    :param filename:
    :param delim: str
    :return:
    """
    insertion_count = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            l_ins, r_ins, l_palindromic, r_palindromic = int(line[header.index('Left.insertion.length')]), \
                int(line[header.index('Right.insertion.length')]), int(line[header.index('Left.palindromic')]), \
                int(line[header.index('Right.palindromic')])
            if (l_ins - l_palindromic) + (r_ins - r_palindromic) != 0:
                insertion_count += 1

    return insertion_count


if __name__ == '__main__':
    files = [
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal.tsv',
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal.tsv'
    ]
    for file in files:
        count = check_insertions(file)
        fn = file.split('\\')[-1]
        print(f'Number of duplicate sequences in file \"{fn}\": {count}')
