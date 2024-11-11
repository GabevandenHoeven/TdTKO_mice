import numpy as np
import csv
import matplotlib.pyplot as plt
import pandas as pd


def read_pgens_from_file(filename, column_index, delim='\t'):
    """

    :param filename:
    :param column_index:
    :param delim:
    :return:
    """
    generation_probabilities = []
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            generation_probabilities.append(float(line[header.index(column_index)]))
    return generation_probabilities


def get_bins(data, n_bins, n_total_sequences):
    """

    :param data:
    :param n_bins:
    :param n_total_sequences:
    :return:
    """
    hist, bins = np.histogram(data, bins=n_bins,)
    print("Bin Edges:", bins)
    print("Histogram Counts:", hist)
    return [(bins[e] + bins[e + 1]) / 2 for e in range(len(hist))], [e / n_total_sequences for e in hist]


def plot_pgens_distribution(plot_title, axislabels, legendlabels, outfile, x_data, frequencies):
    """

    :param plot_title:
    :param axislabels:
    :param legendlabels:
    :param outfile:
    :param x_data:
    :param frequencies:
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

    plt.plot(x_data[0], frequencies[0], label=legendlabels[0], color=okabe_ito_colours[4])
    plt.plot(x_data[1], frequencies[1], label=legendlabels[1], color=okabe_ito_colours[5])
    plt.legend()
    plt.savefig(outfile)
    plt.close()
    return


def filter_zeroes(data, lower_bound=None):
    """

    :param data:
    :param lower_bound:
    :return:
    """
    df_pgens = pd.DataFrame(data)
    df_pgens.columns = ['Pgen']
    print("Old shape:", df_pgens.shape)
    if lower_bound is not None:
        filtered_array = np.where(df_pgens <= lower_bound)[0]
        df_pgens.drop(index=filtered_array, inplace=True)
    print("new shape:", df_pgens.shape)
    zeroes = np.where(df_pgens == 0.0)[0]
    df_pgens.drop(index=zeroes, inplace=True)
    print("Zeroes filtered shape:", df_pgens.shape)
    return df_pgens['Pgen'].tolist()


def distribute_pgens(datafile1, datafile2, column_name):
    x_ = []
    y_ = []
    num_bins = 25

    treg_pgens = read_pgens_from_file(datafile1, column_name)
    naive_pgens = read_pgens_from_file(datafile2, column_name)
    treg_total = len(treg_pgens)
    naive_total = len(naive_pgens)

    treg_pgens = filter_zeroes(treg_pgens, 1.0e-30)
    naive_pgens = filter_zeroes(naive_pgens, 1.0e-30)
    treg_pgens = np.log10(treg_pgens)
    pgens, freq = get_bins(treg_pgens, num_bins, treg_total)
    x_.append(pgens)
    y_.append(freq)

    naive_pgens = np.log10(naive_pgens)
    pgens, freq = get_bins(naive_pgens, num_bins, naive_total)
    x_.append(pgens)
    y_.append(freq)

    print(x_)
    print(y_)
    return x_, y_


if __name__ == '__main__':
    tregs_datafile = '..\\data_files\\Tregs_data\\control\\CD4_Treg\\filtered_data\\filtered_data_Treg_with_pgens.tsv'
    tnaive_datafile = '..\\data_files\\Tregs_data\\control\\CD4_Naive\\filtered_data\\filtered_data_TNaive_with_pgens.tsv'

    title = 'Distribution of Nucleotide Generation Probabilities'
    x_labels = ('nt Pgen', 'Frequency')
    l_labels = ('Tregs', 'Naive')
    out_filename = 'img/pgen_distribution_nt.png'
    x, y = distribute_pgens(tregs_datafile, tnaive_datafile, column_name='nucleotide_pgen')
    plot_pgens_distribution(title, x_labels, l_labels, out_filename, x, y)

    title = 'Distribution of Amino Acid Generation Probabilities'
    x_labels = ('aa Pgen', 'Frequency')
    l_labels = ('Tregs', 'Naive')
    out_filename = 'img/pgen_distribution_aa.png'
    x, y = distribute_pgens(tregs_datafile, tnaive_datafile, column_name='amino_acid_pgen')
    plot_pgens_distribution(title, x_labels, l_labels, out_filename, x, y)
