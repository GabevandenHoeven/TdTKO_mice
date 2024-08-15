import csv


def check_duplicate_cdr3(filename: str, delim='\t'):
    """Checks if the file contains unique CDR3 nucleotide sequences. They are considered unique if the VJ combination
    together with the nucleotide sequence has not been found in the file.

    :param filename:
    :param delim: str
    :return:
    """
    duplicate_count = 0
    unique = {}
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            mouse, v_gene, j_gene, cdr3 = line[header.index('Mouse')], line[header.index('V.gene')], line[header.index('J.gene')], \
                line[header.index('Junction.nucleotide.sequence')]
            try:
                if unique[mouse + v_gene + cdr3 + j_gene]:
                    duplicate_count += 1
            except KeyError:
                unique.update({mouse + v_gene + cdr3 + j_gene: True})
    return duplicate_count


if __name__ == '__main__':
    files = [
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv',
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]
    for file in files:
        count = check_duplicate_cdr3(file)
        fn = file.split('\\')[-1]
        print(f'Number of duplicate sequences in file \"{fn}\": {count}')
