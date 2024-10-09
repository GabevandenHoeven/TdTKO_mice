import numpy as np
from scipy.stats import skew, skewtest, kurtosis, kurtosistest, ttest_ind
from plots import plot_d_lengths
from utils import get_unique_sequences_from_file
from statistics import stdev


def get_d_lengths(datafile):
    """
    :param datafile:
    :return:
    """

    d_lengths_counts = {}
    total = 0
    header = datafile[0]

    for line in datafile[1:]:
        total += 1
        try:
            d_lengths_counts[int(line[header.index('D.length.used')])] += 1
        except KeyError:
            d_lengths_counts.update({int(line[header.index('D.length.used')]): 1})
    return d_lengths_counts, total


def get_d_lengths_per_mouse(datafile):
    """Returns a hashmap of lengths of D segments per mouse.
    :return:
    """
    mice = {
        # mouse : [{D lengths : occurrences}, number of sequences for mouse]
    }
    total = 0
    header = datafile[0]
    for line in datafile[1:]:
        total += 1
        mouse, d = line[header.index('Mouse')], int(line[header.index('D.length.used')])
        if mouse not in mice.keys():
            mice.update({mouse: [{}, 0]})
        if d not in mice[mouse][0].keys():
            mice[mouse][0].update({d: 0})
        mice[mouse][0][d] += 1
        mice[mouse][1] += 1

    return mice, total


if __name__ == '__main__':
    x = []
    y = []
    means = []
    data = []
    files = ['..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
             '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
             ]
    for file in files:
        filtered_data = get_unique_sequences_from_file(file)
        d_lengths, total_lines = get_d_lengths(filtered_data)
        temp_list = [e * d_lengths[e] for e in sorted(d_lengths.keys())]
        means.append(sum(temp_list) / total_lines)
        data_list = []
        for i in range(len(temp_list)):
            data_list.extend([i] * d_lengths[i])
        data.append(data_list)
        print(total_lines, means[-1])

        x.append([e for e in sorted(d_lengths.keys())])
        y.append([d_lengths[d] / total_lines * 100 for d in sorted(d_lengths.keys())])

    st_devs = []
    for i in data:
        # The arguments below are default, though filled in to show the options
        arr = np.array(i)
        skewness, skew_test, kurt, kurt_test = skew(arr), skewtest(arr, alternative='two-sided'), \
            kurtosis(arr, fisher=True), kurtosistest(arr, alternative='two-sided')
        std = stdev(i)
        st_devs.append(std)
        print(skewness, skew_test, kurt, kurt_test, std)
    t_test = ttest_ind(a=np.array(data[0]), b=np.array(data[1]), axis=0, equal_var=True, alternative='two-sided')
    print(f'T test mean TdTKO is different from WT: {t_test}')

    out = f'..\\img\\Distribution_of_D_lengths.png'
    title = f'Distribution of D lengths'
    plot_d_lengths(x, y, means, st_devs, ['TdTKO', 'WT'], title, out)



