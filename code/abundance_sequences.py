import csv


def get_abundance(data, exp: str, max_incidence: int):
    """This function reads a file with CDR3 sequences of mice annotated with V and J identifiers,
    (inferred) D sequence, as well as mouse identifier and (inferred) insertions. It stores the incidence of the
    sequences, then sorts them by incidence, and finally calculates the fraction of sequences per incidence group that
    match a given expression.

    :param data: str - The lines of the file which were first filtered to only contain unique sequences.
    :param exp: str - An expression to check for a specific D length. Write in format 'seq[3] == 0'
    :param max_incidence: int - The number of mice present in the file. If a sequence is present in all mice, this is
    the maximum incidence.
    :returns:
    :return fractions: list - A list with the fractions of sequences that match an expression per incidence.
    :return incidences: list - A list with lists of sequence information per incidence
    """
    sequences = {
        # V_identifier + sequence J_identifier: [count, sequence, [mice], d_length, VJ distance, insertions]
    }
    header = data[0]
    dup = {}
    for line in data[1:]:
        mouse, sequence, d_length, v, j, vj_dis, ins = line[header.index('Mouse')], \
            line[header.index('Junction.nucleotide.sequence')], int(line[header.index('D.length.used')]), \
            line[header.index('V.gene')], line[header.index('J.gene')], \
            int(line[header.index('V.J.distance')]), \
            (int(line[header.index('Left.insertion.length')]) - int(line[header.index('Left.palindromic')]) +
             int(line[header.index('Right.insertion.length')]) - int(line[header.index('Right.palindromic')]))
        try:
            if mouse in sequences[v + sequence + j][2]:
                if 'gen' not in mouse:
                    raise ValueError('This sequence from this mouse is already in the dictionary')
                try:
                    dup[v + sequence + j] += 1
                except KeyError:
                    dup.update({v + sequence + j: 1})
                print(f'Sequence: {v + sequence + j} for mouse {mouse} occurs multiple times.\n'
                      f'Duplicates: {dup[v + sequence + j]}')
                continue
            elif d_length != sequences[v + sequence + j][3]:
                raise ValueError('The same sequence has a different D length')
            elif vj_dis != sequences[v + sequence + j][4]:
                raise ValueError('The same sequence has a different VJ distance')
            elif ins != sequences[v + sequence + j][5]:
                raise ValueError('The same sequence has a different insertion length')
            else:
                sequences[v + sequence + j][0] += 1
                sequences[v + sequence + j][2].append(mouse)
        except KeyError:
            sequences.update({v + sequence + j: [1, sequence, [mouse], d_length, vj_dis, ins, v, j]})

    fractions = []
    incidences = []
    for i in range(1, max_incidence + 1):
        # i is the incidence or number of mice the sequence is found in.
        tmp_list = [c for c in sequences.values() if c[0] == i]
        # Make a list per group of incidence
        incidences.append(tmp_list)

        n_seq_per_incidence = sum(1 for seq in tmp_list if eval(exp))
        # How many sequences in each incidence group match the expression
        n_rows = len(tmp_list)
        try:
            fractions.append(n_seq_per_incidence / n_rows * 100)
        except ZeroDivisionError:
            fractions.append(0)
    return fractions, incidences
