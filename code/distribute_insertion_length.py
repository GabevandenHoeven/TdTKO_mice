from plots import plot_dist_insertion_length
from utils import get_unique_sequences_from_file


def get_insertion_counts(file):
    """Returns a hashmap of the length of insertions and how often they appear together with the mean length
    :param file: nested list - A nested list of a datafile. Starts with a header.
    :return:
    """
    total_seq = 0
    insertions = []
    counts = {}
    ins = None
    header = file[0]
    for row in file[1:]:
        total_seq += 1
        try:
            ins = (int(row[header.index('Left.insertion.length')]) - int(row[header.index('Left.palindromic')])) + \
                  (int(row[header.index('Right.insertion.length')]) - int(row[header.index('Right.palindromic')]))
            insertions.append(ins)
            counts[ins] += 1
        except KeyError:
            counts.update({ins: 1})
    mean = sum(insertions) / total_seq
    return counts, total_seq, mean


if __name__ == '__main__':
    x_points = []
    y_points = []
    labels = ['TdTKO', 'WT']
    averages = []
    file_list = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv',
    ]

    for filename in file_list:
        unique_sequences = get_unique_sequences_from_file(filename)
        insertion_counts, total, average = get_insertion_counts(unique_sequences)
        x = sorted(insertion_counts.keys())
        x_points.append(x)
        y_points.append([insertion_counts[point] / total * 100 for point in x])
        averages.append(average)
        print(f'The average is: {average}')

    outfile = f'..\\img\\Insertion_length_distribution.png'
    plot_dist_insertion_length(x_points, y_points, labels, averages, outfile)
