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

    x_points = [0, 5, 10, 15, 20, 25, 50]

    y_points = [47.94, 21.18, 13.18, 8.90, 6.40, 4.84, 1.89]
    plt.plot(x_points, y_points, label='TdTKO')

    # y_points = [94.29, 39.41, 20.10, 11.39, 6.71, 4.09, 0.54]
    # plt.plot(x_points, y_points, label='mandl-07')

    y_points = [94.11, 36.93, 18.45, 10.68, 6.61, 4.28, 0.78]
    plt.plot(x_points, y_points, label='Normal')

    plt.title('Fraction of sequences with insertions per threshold of supporting reads')
    plt.xlabel('Supporting reads')
    plt.ylabel('Fraction of seq with insertions')
    plt.legend()
    plt.savefig('Normal_TdTKO_fraction_insertions_per_reads.png')
    plt.close()


def plot_sequence_length_read_count(x_points, y_points, label):
    plt.figure()

    plt.scatter(x_points, y_points, label=label, s=3)
    plt.title('Read count per Junction sequence length')
    plt.xlabel('Junction sequence length')
    plt.ylabel('Read count')
    plt.legend()
    plt.savefig(f'Read_count_per_sequence_length_{label}.png')
    plt.close()


if __name__ == '__main__':
    # plot_supporting_reads()
    plot_fractions_ins()