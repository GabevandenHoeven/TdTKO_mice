from utils import get_unique_sequences_from_file
from plots import plot_distribute_pgens


def get_pgens_distribution_data(data):
    header = data[0]
    pgens = {}
    count = 0
    for sequence in data[1:]:
        count += 1
        try:
            pgens[float(sequence[header.index('Generation.probability')])] += 1
        except KeyError:
            pgens.update({float(sequence[header.index('Generation.probability')]): 1})
    return pgens, count


if __name__ == '__main__':

    file = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2 (1).tsv'
    filtered_data = get_unique_sequences_from_file(file)
    generation_probabilities, line_count = get_pgens_distribution_data(filtered_data)
    # TODO this doesnt quite work, think of a different way to show this plot
    #  maybe preset some probability ranges and place probabilities in clusters to get bigger groups

    # plot_distribute_pgens()
