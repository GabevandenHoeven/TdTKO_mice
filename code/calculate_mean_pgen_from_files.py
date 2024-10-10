import csv


def calculate_mean_pgen_from_file(filename, delim='\t'):
    """Calculates the mean pgen from a datafile.

    :param filename: str - The file path of the datafile
    :param delim: the delimiter for the datafile. Default is a tab '\t'
    :return:
    """
    mean_pgen = 0
    number_of_lines = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        for line in reader:
            number_of_lines += 1
            mean_pgen += float(line[1])
        mean_pgen = mean_pgen / number_of_lines
    return mean_pgen


if __name__ == '__main__':
    files = [
        'pgen-out_high_incidence_filtered_data_exp_WT_v2 (1).tsv',
        'pgen-out_low_incidence_filtered_data_exp_WT_v2 (1).tsv',
    ]
    for file in files:
        mean = calculate_mean_pgen_from_file(file)
        print(mean)
