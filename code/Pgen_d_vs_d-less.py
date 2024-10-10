import numpy
from utils import get_unique_sequences_from_file


def get_pgens_for_d_and_d_less_sequences(data):
    pgens_with_d = []
    pgens_without_d = []
    header = data[0]
    for sequence in data[1:]:
        d_length, pgen = int(sequence[header.index('D.length.used')]), \
            float(sequence[header.index('Generation.probability')])
        if d_length == 0:
            pgens_without_d.append(pgen)
        else:
            pgens_with_d.append(pgen)
    return pgens_with_d, pgens_without_d


if __name__ == '__main__':
    file = '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2 (1).tsv'
    filtered_data = get_unique_sequences_from_file(file)
    pgens_d, pgens_no_d = get_pgens_for_d_and_d_less_sequences(filtered_data)
    print('with D segment')
    print(len(pgens_d))
    print('mean:', sum(pgens_d) / len(pgens_d))
    print('median:', numpy.median(pgens_d))
    print('without D segment')
    print(len(pgens_no_d))
    print('mean:', sum(pgens_no_d) / len(pgens_no_d))
    print('median:', numpy.median(pgens_no_d))
    print((len(pgens_no_d) / (len(pgens_d) + len(pgens_no_d))) * 100, '% of sequences')
