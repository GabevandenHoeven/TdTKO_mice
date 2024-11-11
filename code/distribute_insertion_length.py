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
    sethna_data = [
        [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
         [0.8979591901958919, 0.06337268576736017, 0.01933404972611622, 0.0064446832420386796, 0.0010740483149649955,
          0, 0.0010740483149649955, 0, 0, 0, 0, 0.0010740483149649955, 0.0010740483149649955, 0, 0, 0]],
        [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
         [0.10096663856656544, 0.2370461581654641, 0.24926957437565775, 0.1589687877449141, 0.10955961511465846,
          0.04726097766226329, 0.023630439662100558, 0.008592878210030847, 0.0032223416210193398, 0.0010740483149649955,
          0.0010740483149649955, 0.0010740483149649955, 0, 0, 0.0010740483149649955, 0]]
    ]
    # This data was reverse engineered using PlotDigitizer from figure 1 in Sethna, Z., et al. (2017).
    # Insights into immune system development and function from mouse T-cell repertoires.
    # Proceedings of the National Academy of Sciences, 114(9), 2253–2258. https://doi.org/10.1073/pnas.1700241114
    for i in range(len(sethna_data)):
        sethna_data[i][1] = [e * 100 for e in sethna_data[i][1]]
        # transforming frequencies to percentages

    plot_dist_insertion_length(x_points, y_points, sethna_data, labels, averages, outfile)
