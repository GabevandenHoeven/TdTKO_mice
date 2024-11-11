import os
import csv


def filter_data_file(dataframe, filename, delim='\t'):
    """

    :param dataframe:
    :param filename:
    :param delim:
    :return:
    """
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        row_count = 0
        for line in reader:
            if line[header.index('frame_type')] == 'In':
                row_count += 1
                new_line = '\t'.join([str(row_count), line[header.index('sample_name')], line[header.index('species')],
                                      line[header.index('locus')], line[header.index('sample_tags')],
                                      line[header.index('cdr3_rearrangement')], line[header.index('cdr3_amino_acid')],
                                      line[header.index('cdr3_length')]]) + '\n'
                dataframe.append(new_line)
    return dataframe


def read_data_files(data_lines, d, out_fn):
    """

    :param data_lines:
    :param out_fn:
    :param d:
    :return:
    """

    for file in os.listdir(d):
        f = os.path.join(d, file)
        if os.path.isfile(f):
            data_lines = filter_data_file(data_lines, f)

    with open(out_fn, 'w') as outfile:
        outfile.writelines(data)
    return


if __name__ == '__main__':
    data = ['column\tsample_name\tspecies\tlocus\tsample_tags\tnucleotide_sequence\tamino_acid_sequence\tcdr3_length\n']

    directory = '..\\data_files\\Tregs_data\\control\\CD4_Treg'
    output_filename = '..\\data_files\\Tregs_data\\control\\CD4_Treg\\filtered_data\\filtered_data_Tregs.tsv'
    read_data_files(data, directory, output_filename)

    directory = '..\\data_files\\Tregs_data\\control\\CD4_Naive'
    output_filename = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\filtered_data_TNaive.tsv'
    read_data_files(data, directory, output_filename)

