import csv
import matplotlib.pyplot as plt


def read_pgens_from_file(filename, max_incidence, columns, delim='\t'):
    incidence_data = {}
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            sample, seq, pgen = line[header.index('sample_name')], line[header.index(columns[0])], \
                line[header.index(columns[1])]
            try:
                incidence_data[seq][1].append(sample)
            except KeyError:
                incidence_data.update({seq: [pgen, [sample]]})
    incidences = []
    for i in range(1, max_incidence + 1):
        # i is the number of humans a sequence is found in
        tmp_list = [float(incidence_data[seq][0]) for seq in incidence_data.keys() if len(incidence_data[seq][1]) == i]
        incidences.append(tmp_list)
    return incidences


def boxplot_pgens_incidence(plot_title, axislabels, legendlabels, outfile, data, pos):
    """

    :param plot_title:
    :param axislabels:
    :param legendlabels:
    :param outfile:
    :param data:
    :param pos:
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
    plt.yscale('log')
    width = [0.4] * len(pos)
    bplot = plt.boxplot(data[0], positions=[e - 0.25 for e in pos], showfliers=False, showmeans=True, meanline=True,
                        widths=width, patch_artist=True, vert=True, manage_ticks=False)
    label_boo = True
    for patch in bplot['boxes']:
        patch.set_facecolor(okabe_ito_colours[4])
        if label_boo:
            patch.set_label(legendlabels[0])
            label_boo = False
    bplot2 = plt.boxplot(data[1], positions=[e + 0.25 for e in pos], showfliers=False, showmeans=True, meanline=True,
                         widths=width, patch_artist=True, vert=True, manage_ticks=False)
    label_boo = True
    for patch in bplot2['boxes']:
        patch.set_facecolor(okabe_ito_colours[5])
        if label_boo:
            patch.set_label(legendlabels[1])
            label_boo = False
    plt.axhline(0, 0.5, 0.5, linestyle='dashed', color='green', label='Mean')
    plt.axhline(0, 0.5, 0.5, color='orange', label='Median')
    plt.legend()
    plt.savefig(outfile)
    plt.close()
    return


if __name__ == '__main__':
    tregs_datafile = '..\\data_files\\Tregs_data\\control\\CD4_Treg\\filtered_data\\filtered_data_Treg_with_pgens.tsv'
    tnaive_datafile = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\filtered_data_TNaive_with_pgens.tsv'
    incidence = 8
    # column_names = ('nucleotide_sequence', 'nucleotide_pgen')
    column_names = ('amino_acid_sequence', 'amino_acid_pgen')
    treg_pgens_per_incidence = read_pgens_from_file(tregs_datafile, incidence, column_names)
    naive_pgens_per_incidence = read_pgens_from_file(tnaive_datafile, incidence, column_names)

    # title = 'Nucleotide Generation Probabilities per Incidence'
    title = 'Amino Acid Generation Probabilities per Incidence'
    # x_labels = ('Incidence', 'nt Pgen')
    x_labels = ('Incidence', 'aa Pgen')
    # out_filename = 'img/nt_pgens_of_tregs_tnaive_boxplot.png'
    out_filename = 'img/aa_pgens_of_tregs_tnaive_boxplot.png'
    l_labels = ('Tregs', 'Naive')
    boxplot_pgens_incidence(title, x_labels, l_labels, out_filename,
                            [treg_pgens_per_incidence, naive_pgens_per_incidence],
                            [i for i in range(1, incidence + 1)])
