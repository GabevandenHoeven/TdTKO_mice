import csv
import matplotlib.pyplot as plt


def read_cdr3_length_from_file(filename, delim='\t'):

    cdr3_lengths = {}
    sequences = {}
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            nt_seq, cdr3_length = line[header.index('nucleotide_sequence')], int(line[header.index('cdr3_length')])
            try:
                if sequences[nt_seq]:
                    pass
            except KeyError:
                sequences.update({nt_seq: True})
                # Only unique sequences. This is slightly flawed because we don't have the segments called by RTCR.
                # This means that some sequences are actually unique but just happen to get the same CDR3 sequence.
                try:
                    cdr3_lengths[cdr3_length] += 1
                except KeyError:
                    cdr3_lengths.update({cdr3_length: 1})
    return cdr3_lengths


def plot_cdr3_lengths(plot_title, axislabels, legendlabels, outfile, data1, data2, means):
    """

    :param plot_title:
    :param axislabels:
    :param legendlabels:
    :param outfile:
    :param data1:
    :param data2:
    :param means:
    :return:
    """

    okabe_ito_colours = ['#000000', '#E69F00', '#56B4E9',
                         '#009E73', '#F0E442', '#0072B2',
                         '#D55E00', '#CC79A7', '#999999']
    # 0: black, 1: light orange, 2: light blue, 3: teal, 4: yellow, 5: dark blue, 6: dark orange, 7: pink, 8: gray
    plt.figure()
    plt.title(plot_title)
    plt.xlabel(axislabels[0])
    plt.ylabel(axislabels[1])

    plt.plot(data1[0], data1[1], label=legendlabels[0], color=okabe_ito_colours[4])
    plt.scatter(data1[0], data1[1], s=10, color=okabe_ito_colours[4])
    plt.plot(data2[0], data2[1], label=legendlabels[1], color=okabe_ito_colours[5])
    plt.scatter(data2[0], data2[1], s=10, color=okabe_ito_colours[5])

    plt.axvline(means[0], linestyle='dashed', color=okabe_ito_colours[4], label='mean')
    plt.axvline(means[1], linestyle='dashed', color=okabe_ito_colours[5], label='mean')
    plt.legend()
    # plt.show()
    plt.savefig(outfile)
    return


if __name__ == '__main__':
    tregs_datafile = '..\\data_files\\Tregs_data\\control\\CD4_Treg\\filtered_data\\filtered_data_Treg_with_pgens.tsv'
    tnaive_datafile = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\filtered_data_TNaive_with_pgens.tsv'

    treg_cdr3_lengths = read_cdr3_length_from_file(tregs_datafile)
    naive_cdr3_lengths = read_cdr3_length_from_file(tnaive_datafile)

    cdr3_data1 = [treg_cdr3_lengths[e] for e in sorted(treg_cdr3_lengths.keys())]
    cdr3_data1_keys = [e for e in sorted(treg_cdr3_lengths.keys())]
    cdr3_data2 = [naive_cdr3_lengths[e] for e in sorted(naive_cdr3_lengths.keys())]
    cdr3_data2_keys = [e for e in sorted(naive_cdr3_lengths.keys())]

    mean = [sum([e * treg_cdr3_lengths[e] for e in treg_cdr3_lengths.keys()]) / sum(cdr3_data1),
            sum([e * naive_cdr3_lengths[e] for e in naive_cdr3_lengths.keys()]) / sum(cdr3_data2)]

    cdr3_data1 = [e / sum(cdr3_data1) for e in cdr3_data1]
    cdr3_data2 = [e / sum(cdr3_data2) for e in cdr3_data2]

    cdr3_data1 = [cdr3_data1_keys, cdr3_data1]
    cdr3_data2 = [cdr3_data2_keys, cdr3_data2]
    title = 'Distribution of nt CDR3 lengths'
    x_labels = ('CDR3 length (nt)', 'Frequency')
    l_labels = ('Tregs', 'Naive')
    out_filename = 'img/cdr3_lengths_of_tregs_tnaive.png'

    plot_cdr3_lengths(title, x_labels, l_labels, out_filename, cdr3_data1, cdr3_data2, mean)
