from plots import plot_fractions_ins
import csv


def get_fraction_per_read_threshold(filename: str, thresholds: list, delim='\t'):
    """

    :param filename:
    :param thresholds:
    :param delim:
    :return:
    """
    y = []
    for threshold in thresholds:
        count_lines_supporting_reads = 0
        count_ins = 0
        with open(filename, 'r') as in_file:
            reader = csv.reader(in_file, delimiter=delim)
            header = next(reader)
            for line in reader:
                if int(line[header.index('Number.of.reads')]) > threshold:
                    count_lines_supporting_reads += 1
                    ins = (int(line[header.index('Left.insertion.length')]) -
                           int(line[header.index('Left.palindromic')])) + \
                          (int(line[header.index('Right.insertion.length')]) -
                           int(line[header.index('Right.palindromic')]))
                    if ins > 0:
                        count_ins += 1
        fraction = count_ins / count_lines_supporting_reads * 100
        y.append(fraction)
    return y


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal.tsv',
    ]
    thres = [0, 1, 2, 5, 10, 15, 20, 25, 50]
    y_points = []
    for file in files:
        fracts = get_fraction_per_read_threshold(file, thres)
        y_points.append(fracts)

    title = '..\\img\\Normal-TdTKO_fraction_insertions_per_reads (1).png'
    plot_fractions_ins(thres, y_points, title)
