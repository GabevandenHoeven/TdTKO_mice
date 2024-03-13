import numpy
from matplotlib import pyplot as plt


def plot_supporting_reads():
    """
    """
    with open('count_reads.txt', 'r') as file:
        header = next(file)
        x_points = []
        y_points = []
        for row in file:
            x_points.append(int(row.split(': ')[0]))
            y_points.append(int(row.split(': ')[1].rstrip('\n')))
        x_ = [1, 3, 5, 8, 10, 15, 20, 25, 50, 75, 100, 150, 200, 250, 300]
        y_ = [y_points[e] for e in x_]

        plt.figure()
        plt.plot(x_, y_)
        plt.title('Number of CDR3 sequences supported by number of reads')
        plt.ylabel('Sequences')
        plt.xlabel('Reads')
        plt.savefig('sequences_per_reads.png')
        plt.close()


def plot_fractions_ins():
    plt.figure()

    # old ----------------------------------------------------------------------------------
    # x_points = [0, 5, 10, 15, 20, 25, 50]
    # y_points = [47.94, 21.18, 13.18, 8.90, 6.40, 4.84, 1.89]
    # plt.plot(x_points, y_points, label='TdTKO')

    # y_points = [94.29, 39.41, 20.10, 11.39, 6.71, 4.09, 0.54]
    # plt.plot(x_points, y_points, label='mandl-07')

    # y_points = [94.11, 36.93, 18.45, 10.68, 6.61, 4.28, 0.78]
    # plt.plot(x_points, y_points, label='Normal')
    # -------------------------------------------------------------------------------------
    x_points = [0, 1, 2, 5, 10, 15, 20, 25, 50]

    y_points = [47.94, 36.92, 35.54, 32.07, 27.96, 25.08, 22.90, 21.25, 16.33]
    plt.plot(x_points, y_points, label='TdTKO')

    y_points = [94.11, 92.65, 91.87, 90.33, 88.60, 86.98, 85.33, 83.66, 74.89]
    plt.plot(x_points, y_points, label='Normal')

    plt.title('Fraction of sequences with insertions per threshold of supporting reads')
    plt.xlabel('Supporting reads')
    plt.ylabel('Fraction of seq with insertions')
    plt.xticks([0, 1, 2, 5, 10, 15, 20, 25, 50])
    plt.legend()
    plt.savefig('C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\Normal_TdTKO_fraction_insertions_per_reads.png')
    plt.close()


def plot_sequence_length_read_count(x_points, y_points, label):
    plt.figure()

    plt.scatter(x_points, y_points, label=label, s=3)
    plt.title('Read count per Junction sequence length')
    plt.xlabel('Junction sequence length (nt)')
    plt.ylabel('Read count')
    plt.legend()
    plt.savefig(f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\Read_count_per_sequence_length_{label}.png')
    plt.close()


def plot_dist_junction_sequence_length(x_points_lists, y_points_lists, labels, averages):
    plt.figure()
    x_points, y_points, label, avg = x_points_lists[0], y_points_lists[0], labels[0], averages[0]
    plt.plot(x_points, y_points, label=label)
    plt.axvline(avg, linestyle='dashed', label='average '+label)

    x_points, y_points, label, avg = x_points_lists[1], y_points_lists[1], labels[1], averages[1]
    plt.plot(x_points, y_points, label=label)
    plt.axvline(avg, linestyle='dashed', color='orange', label='average ' + label)

    plt.title('Distribution of junction sequence lengths')
    plt.xlabel('Junction sequence length (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(
        f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\Junction_length_distribution_{"-".join(labels)}')
    plt.close()


def plot_dist_insertion_length(x_points_lists, y_points_lists, labels):
    plt.figure()
    for e in range(len(x_points_lists)):
        x_points, y_points, label = x_points_lists[e], y_points_lists[e], labels[e]
        plt.plot(x_points, y_points, label=label)
    plt.title('Distribution of insertion lengths')
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.xlabel('Insertion length (nt)')
    plt.ylabel('Percentage of sequences with an insertion (%)')
    plt.legend()
    plt.savefig(
        f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\Insertion_length_distribution_{"-".join(labels)}')
    plt.close()


def plot_vj_distance_reads(x_points_lists, y_points_lists, labels):
    plt.figure()
    for e in range(len(x_points_lists)):
        x_points, y_points, label = x_points_lists[e], y_points_lists[e], labels[e]
        plt.scatter(x_points, y_points, label=label, s=3)
    plt.title('Distribution of VJ distances')
    plt.xlabel('VJ distance (nuc)')
    plt.ylabel('Number of reads')
    plt.legend()
    plt.savefig(
        f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\VJ_distance_distribution_{"-".join(labels)}')
    plt.close()


def plot_vj_distance_perc(x_points_lists, y_points_lists, avg, labels):
    plt.figure()

    x_points, y_points, label, average = x_points_lists[0], y_points_lists[0], labels[0], avg[0]
    plt.plot(x_points, y_points, label=label)
    plt.axvline(average, linestyle='dashed', label='average '+label)

    x_points, y_points, label, average = x_points_lists[1], y_points_lists[1], labels[1], avg[1]
    plt.plot(x_points, y_points, label=label)
    plt.axvline(average, linestyle='dashed', color='orange', label='average ' + label)

    plt.title('Distribution of VJ distances')
    plt.xlabel('VJ distance (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(
        f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\VJ_distance_distribution_seq_fraction_{"-".join(labels)}')
    plt.close()


def plot_deletions(x_points_list, y_points_list, labels, title, outfile):
    plt.figure()
    x = x_points_list[0]
    x.extend(x_points_list[1])
    x = numpy.asarray(list(set(x)))
    width = 0.40
    y_points, label = y_points_list[0], labels[0]
    plt.bar(x-0.2, y_points, width=width, label=label)
    y_points, label = y_points_list[1], labels[1]
    plt.bar(x+0.2, y_points, width=width, label=label)
    plt.xticks(numpy.arange(24))
    plt.title(title)
    plt.xlabel('Deletions (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(outfile)
    plt.close()


if __name__ == '__main__':
    # plot_supporting_reads()
    plot_fractions_ins()
