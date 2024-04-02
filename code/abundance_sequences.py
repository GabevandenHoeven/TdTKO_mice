import csv
from plots import plot_line_and_scatter_per_incidence


def get_abundance(filename: str, exp: str):
    """

    :param filename: str - The name of the file with sequences
    :param exp: str - An expression to check for a specific D length
    :return:
    """
    sequences = {
        # sequence: [count, [d_lengths], [mice]]
    }
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)
        for line in reader:
            mouse, sequence, d_length, v, j, pheno, vj_dis, ins = line[header.index('mouse')], \
                line[header.index('Junction.nucleotide.sequence')], int(line[header.index('D.length.used')]), \
                line[header.index('V.gene')], line[header.index('J.gene')], line[header.index('phenotype')], \
                int(line[header.index('V.J.distance')]), int(line[header.index('insertion.length')])
            try:
                if eval(exp) and pheno == 'CD4+':
                    sequences[v + sequence + j][0] += 1
                    sequences[v + sequence + j][1].append(d_length)
                    sequences[v + sequence + j][2].append(mouse)
                    sequences[v + sequence + j][3].append(vj_dis)
                    sequences[v + sequence + j][4].append(ins)
            except KeyError:
                sequences.update({v + sequence + j: [1, [d_length], [mouse], [vj_dis], [ins]]})
    with open(filename, 'r') as file:
        next(file)
        n_rows = len(file.readlines())
    return sequences, n_rows


if __name__ == '__main__':
    files = [
        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv',
        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
    ]
    mice_per_file = [13, 10]
    x = []
    y = []

    for file in files:
        x.append([i for i in range(1, mice_per_file[files.index(file)] + 1)])
        seqs, n_rows = get_abundance(file, 'd_length >= 0')
        mean_d_lengths = []
        for i in range(1, mice_per_file[files.index(file)] + 1):
            n_seq_per_incidence = sum(1 for c in seqs.values() if c[0] == i)
            # i is the incidence or number of mice the sequence is found in.
            # Here you get the number of sequences matching the expression per incidence.
            try:
                mean_d = sum(c[1][0] for c in seqs.values() if c[0] == i) / n_seq_per_incidence
                mean_d_lengths.append(mean_d)
            except ZeroDivisionError:
                mean_d = 0

        y.append(mean_d_lengths)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Incidence', 'Inferred D-segment length (nt)'],
                                           'Mean inferred D-length',
                                           'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                           '\\Mean_inferred_D_length_per_incidence.png')
    x = []
    y = []
    # TODO EVEN WITH THE VJ INCLUDED IN THE SEQUENCE KEY, WE GET MULTIPLE OCCURRENCES IN THE FILE FROM THE SAME MICE
    #  INTERESTINGLY IT'S USUALLY * 3.
    #  UPDATE: IT'S IN THE DATA, SEQUENCES ARE COMING FROM CD4+ CD5HI AND CD5LO. PROBABLY I SHOULD FILTER FOR CD4+.
    #  I FILTERED FOR CD4+, I GOT 2588 IN SEQS BEFORE, NOW I GET 1648
    #  it seems to go well, edit the plots so that normal only runs till 10 and we have scatter plus line with blue
    #  and orange
    for file in files:
        x.append([i for i in range(1, mice_per_file[files.index(file)] + 1)])
        seqs, n_rows = get_abundance(file, 'd_length == 0')
        incidences = []
        for i in range(1, mice_per_file[files.index(file)] + 1):
            n_seq_per_incidence = sum(1 for c in seqs.values() if c[0] == i)
            incidences.append(n_seq_per_incidence / n_rows * 100)
        y.append(incidences)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Incidence', 'Percentage of sequences (%)'],
                                           'Inferred D-length = 0 nt',
                                           'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                           '\\Fraction_no_D_per_incidence.png')
    x = []
    y = []
    for file in files:
        x.append([i for i in range(1, mice_per_file[files.index(file)] + 1)])
        seqs, n_rows = get_abundance(file, 'd_length <= 2')
        incidences = []
        for i in range(1, mice_per_file[files.index(file)] + 1):
            n_seq_per_incidence = sum(1 for c in seqs.values() if c[0] == i)
            incidences.append(n_seq_per_incidence / n_rows * 100)
        y.append(incidences)
    plot_line_and_scatter_per_incidence(x, y, ['TdTKO', 'Normal'], ['Incidence', 'Percentage of sequences (%)'],
                                           'Inferred D-length <= 2 nt',
                                           'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img'
                                           '\\Fraction_max_2nt_D_per_incidence.png')
