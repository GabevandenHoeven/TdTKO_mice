from utils import get_unique_sequences_from_file
from plots import plot_distribute_pgens


def get_d_length_and_pgen(data):
    """This function takes a nested list of a datafile that contains D segment lengths and the generation probability
    of their sequences. Sorts this data into bins of Pgens per D segment length and returns a hashmap.

    :param data: nested list - A nested list of a datafile. Starts with a header.
    :return:
    """
    res = {}
    header = data[0]
    for line in data[1:]:
        d_length, pgen = int(line[header.index('D.length.used')]), float(line[header.index('Generation.probability')])
        try:
            res[d_length].append(pgen)
        except KeyError:
            res.update({d_length: [pgen]})
    return res


if __name__ == '__main__':
    file = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2 (1).tsv'
    filtered_data = get_unique_sequences_from_file(file)
    results = get_d_length_and_pgen(filtered_data)
    x = [d for d in sorted(results.keys())]
    mean_pgens = [sum(results[d]) / len(results[d]) for d in sorted(results.keys())]
    sequences_per_length = [len(results[d]) for d in sorted(results.keys())]
    pgen_data = [results[d] for d in sorted(results.keys())]
    print(mean_pgens)
    print(sequences_per_length)
    out_file = '..\\img\\Pgen_per_d_length_boxplot.png'
    positions = [i for i in range(0, len(pgen_data))]
    plot_distribute_pgens(pgen_data, positions, 'Generation probabilities per D segment length',
                          ['D segment length (nt)', 'Pgen'], out_file)
