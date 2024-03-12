import pandas
from plots import plot_vj_distance_reads, plot_vj_distance_occ


def calculate_average_v_j_distance(fn, delim):
    """

    :param fn:
    :param delim:
    :return:
    """
    df = pandas.read_csv(fn, sep=delim, header=0)
    sum_ = df['V.J.distance'].sum()
    reads = list(df['Number.of.reads'])
    distances = list(df['V.J.distance'])
    distance_occ = list(df['V.J.distance'].value_counts().sort_index())
    unique_distances = sorted(list(df['V.J.distance'].unique()))
    n_rows = len(df)
    avg = sum_ / n_rows
    print(f'The average VJ distance is: {sum_} / {n_rows} = {avg}')
    return distances, reads, unique_distances, distance_occ


if __name__ == '__main__':
    vj_distances = []
    read_counts = []
    unique_distances = []
    distance_occurrences = []
    labels = []
    filenames = [
        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv',
        'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
                ]
    for file in filenames:
        d, r, unique_d, d_occ = calculate_average_v_j_distance(file, '\t')
        vj_distances.append(d)
        read_counts.append(r)
        unique_distances.append(unique_d)
        d_occ = [occ * 100 / sum(d_occ) for occ in d_occ]
        distance_occurrences.append(d_occ)
        labels.append(file.split('_')[-1].rstrip('.tsv'))
    # plot_vj_distance_reads(vj_distances, read_counts, labels)
    plot_vj_distance_occ(unique_distances, distance_occurrences, labels)
