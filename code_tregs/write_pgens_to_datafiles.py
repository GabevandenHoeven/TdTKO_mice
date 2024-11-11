import csv


def read_pgens_from_output_file(filename: str, data_frame: list):

    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for line in reader:
            data_frame.append((line[1], line[3]))
    return data_frame


def write_pgens_to_datafile(filename, pgen_data, delim='\t'):

    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        new_header = [e for e in header]
        new_header.append('nucleotide_pgen')
        new_header.append('amino_acid_pgen')
        outfile = ['\t'.join(new_header) + '\n']
        for line in reader:
            try:
                pgens = pgen_data[reader.line_num - 2]
            except IndexError:
                break
            line.append(pgens[0])
            line.append(pgens[1])
            outfile.append('\t'.join(line) + '\n')
    filename = filename.rstrip('.tsv') + '_with_pgens.tsv'
    with open(filename, 'w') as out_file:
        out_file.writelines(outfile)
    return


if __name__ == '__main__':
    tregs_pgens_file = '..\\data_files\\Tregs_data\\control\\CD4_Treg\\filtered_data\\pgen_output.tsv'
    tnaive_pgens_file1 = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\pgen_output.tsv'
    tnaive_pgens_file2 = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\pgen_output (1).tsv'
    tnaive_pgens_file3 = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\pgen_output (2).tsv'
    tregs_datafile = '..\\data_files\\Tregs_data\\control\\CD4_Treg\\filtered_data\\filtered_data_Tregs.tsv'
    tnaive_datafile = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\filtered_data_TNaive.tsv'
    pgens = []
    pgens = read_pgens_from_output_file(tregs_pgens_file, pgens)
    write_pgens_to_datafile(tregs_datafile, pgens)

    pgens = []
    pgens = read_pgens_from_output_file(tnaive_pgens_file1, pgens)
    pgens = read_pgens_from_output_file(tnaive_pgens_file2, pgens)
    pgens = read_pgens_from_output_file(tnaive_pgens_file3, pgens)
    write_pgens_to_datafile(tnaive_datafile, pgens)
