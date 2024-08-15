from plots import plot_vj_distance_perc
from utils import get_unique_sequences_from_file


def calculate_average_v_j_distance(lines):
    """

    :param lines:
    :return:
    """
    header = lines[0]
    distances = {}
    total = 0
    for line in lines[1:]:
        total += 1
        try:
            distances[int(line[header.index('V.J.distance')])] += 1
        except KeyError:
            distances.update({int(line[header.index('V.J.distance')]): 1})
    return distances, total


if __name__ == '__main__':
    x = []
    y = []
    labels = ['TdTKO', 'Normal']
    averages = []
    filenames = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
        # '..\\data_files\\TdTKo\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
                ]
    for file in filenames:
        unique_lines = get_unique_sequences_from_file(file)
        vj_distances, total_lines = calculate_average_v_j_distance(unique_lines)
        x.append([e for e in sorted(vj_distances.keys())])
        y_ = [vj_distances[e] / total_lines * 100 for e in sorted(vj_distances.keys())]
        y.append(y_)
        average = sum([e * vj_distances[e] for e in vj_distances.keys()]) / total_lines
        averages.append(average)
    out_fn = f'..\\img\\unique_seq_img\\VJ_distance_distribution_seq_fraction_{"-".join(labels)}'
    # out_fn = f'..\\img\\unique_seq_img\\Generated_VJ_distance_distribution_seq_fraction_{"-".join(labels)}'
    # plot_vj_distance_perc(x, y, averages, labels, out_fn)
    print(averages)
