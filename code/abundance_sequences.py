import csv
from plots import plot_line_and_scatter_per_incidence


def get_abundance(filename: str):
    """
    :param filename: str - The name of the file with sequences

    :return:
    """
    sequences = {
        # sequence: [count, [d_lengths], [mice], [VJ distances], [insertions]]
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
                if mouse in sequences[v + sequence + j][2]:
                    print('This sequence from this mouse is already in the dictionary')
                else:
                    sequences[v + sequence + j][0] += 1
                    sequences[v + sequence + j][1].append(d_length)
                    sequences[v + sequence + j][2].append(mouse)
                    sequences[v + sequence + j][3].append(vj_dis)
                    sequences[v + sequence + j][4].append(ins)
            except KeyError:
                sequences.update({v + sequence + j: [1, [d_length], [mouse], [vj_dis], [ins]]})


    return

