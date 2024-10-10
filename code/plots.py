import numpy
from matplotlib import pyplot as plt
from utils import calculate_confidence_intervals
import scipy.stats as stats


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


def plot_fractions_ins(x_points, y_points, title, labels):
    plt.figure()

    # old ----------------------------------------------------------------------------------
    # x_points = [0, 5, 10, 15, 20, 25, 50]
    # y_points = [47.94, 21.18, 13.18, 8.90, 6.40, 4.84, 1.89]
    # plt.plot(x_points, y_points, label='TdTKO')

    # y_points = [94.29, 39.41, 20.10, 11.39, 6.71, 4.09, 0.54]
    # plt.plot(x_points, y_points, label='mandl-07')

    # y_points = [94.11, 36.93, 18.45, 10.68, 6.61, 4.28, 0.78]
    # plt.plot(x_points, y_points, label='WT')
    # before p-nt -------------------------------------------------------------------------
    # x_points = [0, 1, 2, 5, 10, 15, 20, 25, 50]
    # y_points = [47.94, 36.92, 35.54, 32.07, 27.96, 25.08, 22.90, 21.25, 16.33] TdTKO
    # y_points = [94.11, 92.65, 91.87, 90.33, 88.60, 86.98, 85.33, 83.66, 74.89] WT
    # -------------------------------------------------------------------------------------
    # after p-nt TdTKO: [13.31, 13.31, 12.41, 10.72, 8.75, 7.41, 6.39, 5.61, 3.75]
    # after p-nt WT: [86.75, 86.75, 85.49, 83.25, 80.69, 78.25, 75.78, 73.24, 60.40]

    plt.plot(x_points, y_points[0], label=labels[0], color='blue')
    plt.scatter(x_points, y_points[0], color='blue', s=10)

    plt.plot(x_points, y_points[1], label=labels[1], color='orange')
    plt.scatter(x_points, y_points[1], color='orange', s=10)

    plt.title('Fraction of sequences with insertions \nper threshold of supporting reads')
    plt.xlabel('Supporting reads')
    plt.ylabel('Percentage of sequences with insertions (%)')
    plt.legend(loc='center right')
    plt.savefig(title)
    plt.close()


def plot_sequence_length_read_count(x_points, y_points, label):
    plt.figure()

    plt.scatter(x_points, y_points, label=label, s=3)
    plt.title('Read count per Junction sequence length')
    plt.xlabel('Junction sequence length (nt)')
    plt.ylabel('Read count')
    plt.legend()
    plt.savefig(f'..\\img\\Read_count_per_sequence_length_{label}.png')
    plt.close()


def plot_dist_junction_sequence_length(x_points_lists, y_points_lists, labels, averages):
    plt.figure()
    x_points, y_points, label, avg = x_points_lists[0], y_points_lists[0], labels[0], averages[0]
    plt.plot(x_points, y_points, label=label, color='blue')
    plt.scatter(x_points, y_points, color='blue', s=10)
    plt.axvline(avg, linestyle='dashed', label='average '+label, color='blue')

    x_points, y_points, label, avg = x_points_lists[1], y_points_lists[1], labels[1], averages[1]
    plt.plot(x_points, y_points, label=label, color='orange')
    plt.scatter(x_points, y_points, color='orange', s=10)
    plt.axvline(avg, linestyle='dashed', color='orange', label='average ' + label)

    plt.title('Distribution of junction sequence lengths')
    plt.xlabel('Junction sequence length (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(
        f'..\\img\\Junction_length_distribution_{"-".join(labels)}')
    plt.close()


def plot_dist_insertion_length(x_points_lists, y_points_lists, labels, averages, title):
    plt.figure()
    x_, y_ = [x for x in x_points_lists[0] if x_points_lists[0].index(x) <= 20], [y for y in y_points_lists[0]
                                                                                  if y_points_lists[0].index(y) <= 20]

    plt.plot(x_, y_, label=labels[0], color='blue')
    plt.scatter(x_, y_, color='blue', s=10)
    plt.axvline(averages[0], linestyle='dashed', color='blue', label='average ' + labels[0])
    x_, y_ = [x for x in x_points_lists[1] if x_points_lists[1].index(x) <= 20], [y for y in y_points_lists[1]
                                                                                  if y_points_lists[1].index(y) <= 20]
    plt.plot(x_, y_, label=labels[1], color='orange')
    plt.scatter(x_, y_, color='orange', s=10)
    plt.axvline(averages[1], linestyle='dashed', color='orange', label='average ' + labels[1])
    plt.xticks(numpy.arange(0, 21))
    plt.title('Distribution of insertion lengths')
    plt.xlabel('Insertion length (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(title)
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
        f'..\\img\\VJ_distance_distribution_{"-".join(labels)}')
    plt.close()


def plot_vj_distance_perc(x_points_lists, y_points_lists, avg, labels, out_file):
    plt.figure()
    x_points = [x for x in x_points_lists[0] if x <= 25]
    y_points, label, average = [y for y in y_points_lists[0] if y_points_lists[0].index(y) <= 25], labels[0], avg[0]
    plt.plot(x_points, y_points, label=label, color='blue')
    plt.scatter(x_points, y_points, color='blue', s=10)
    plt.axvline(average, linestyle='dashed', label='average '+label, color='blue')

    x_points = [x for x in x_points_lists[1] if x <= 25]
    y_points, label, average = [y for y in y_points_lists[1] if y_points_lists[1].index(y) <= 25], labels[1], avg[1]
    plt.plot(x_points, y_points, label=label, color='orange')
    plt.scatter(x_points, y_points, color='orange', s=10)
    plt.axvline(average, linestyle='dashed', color='orange', label='average ' + label)

    plt.title('Distribution of VJ distances')
    plt.xlabel('VJ distance (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(out_file)
    plt.close()


def plot_deletions(x_points_list, y_points_list, labels, ticks, title, outfile):
    plt.figure()
    x = x_points_list[0]
    x.extend(x_points_list[1])
    x = set(x)
    x = numpy.asarray([e for e in x if e <= 13])
    width = 0.40
    y_points, label = [y for y in y_points_list[0] if y_points_list[0].index(y) <= 13], labels[0]
    plt.bar(x-0.2, y_points, width=width, label=label, color='blue')
    y_points, label = [y for y in y_points_list[1] if y_points_list[1].index(y) <= 13], labels[1]
    plt.bar(x+0.2, y_points, width=width, label=label, color='orange')
    plt.xticks(ticks)
    plt.title(title)
    plt.xlabel('Deletions (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(outfile)
    plt.close()


def plot_d_deletions(x_points_list, y_points_list, labels, title, outfile):
    plt.figure()
    x = x_points_list[0]
    x.extend(x_points_list[1])
    x = numpy.asarray(list(set(x)))
    width = 0.40
    y_points, label = y_points_list[0], labels[0]
    plt.bar(x-0.2, y_points, width=width, label=label, color='blue')
    y_points, label = y_points_list[1], labels[1]
    plt.bar(x+0.2, y_points, width=width, label=label, color='orange')
    ticks = list(x)
    tick_labels = [e for e in ticks]
    tick_labels[-1] = 'Unidentified'
    plt.xticks(ticks=ticks, labels=tick_labels)
    plt.title(title)
    plt.xlabel('Deletions (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(outfile)
    plt.close()


def plot_d_lengths(x_points_list, y_points_list, means_, st_devs_, labels, title, outfile):
    plt.figure()
    x = x_points_list[0]
    x.extend(x_points_list[1])
    x = numpy.asarray(list(set(x)))
    width = 0.40
    y_points, label = y_points_list[0], labels[0]
    plt.bar(x-0.2, y_points, width=width, label=label, color='blue')
    plt.axvline(means_[0], linestyle='dashed', color='blue', label='average ' + label)

    n = 100
    norm_ = numpy.linspace(start=means_[0] - 3 * st_devs_[0], stop=means_[0] + 3 * st_devs_[0], num=n)
    pdf = stats.norm.pdf(norm_, means_[0], st_devs_[0])
    plt.plot(norm_, n * pdf, color='blue')

    y_points, label = y_points_list[1], labels[1]
    plt.bar(x+0.2, y_points, width=width, label=label, color='orange')
    plt.axvline(means_[1], linestyle='dashed', color='orange', label='average ' + label)

    norm_ = numpy.linspace(start=means_[1] - 3 * st_devs_[1], stop=means_[1] + 3 * st_devs_[1], num=n)
    pdf = stats.norm.pdf(norm_, means_[1], st_devs_[1])
    plt.plot(norm_, n * pdf, color='orange')

    plt.xticks(numpy.arange(15))
    plt.title(title)
    plt.xlabel('D length (nt)')
    plt.ylabel('Percentage of sequences (%)')
    plt.legend()
    plt.savefig(outfile)
    plt.close()


def plot_line_and_scatter_per_incidence(x_points_list, y_points_list, labels, plot_labels, title, out):
    plt.figure()
    plt.plot([n / len(x_points_list[0]) for n in x_points_list[0]], y_points_list[0], label=labels[0], color='blue')
    plt.scatter([n / len(x_points_list[0]) for n in x_points_list[0]], y_points_list[0], color='blue', s=10)
    plt.plot([n / len(x_points_list[1]) for n in x_points_list[1]], y_points_list[1], label=labels[1], color='orange')
    plt.scatter([n / len(x_points_list[1]) for n in x_points_list[1]], y_points_list[1], color='orange', s=10)
    plt.xticks(numpy.arange(0, 1.1, 0.1))
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.title(title)
    plt.legend()
    plt.savefig(out)
    plt.close()


def plot_boxplot_per_incidence(data, labels, positions, plot_labels, ticks, tick_labels, title, out):
    plt.figure(figsize=(10, 7))
    bplot1 = plt.boxplot(data[0], showfliers=False, patch_artist=True, showmeans=True, positions=positions[0])
    bplot2 = plt.boxplot(data[1], showfliers=False, patch_artist=True, showmeans=True, positions=positions[1])
    colors = ['blue', 'orange']
    for bplot in (bplot1, bplot2):
        color = colors[(bplot1, bplot2).index(bplot)]
        for patch in bplot['boxes']:
            patch.set_facecolor(color)
        for median in bplot['medians']:
            if color == 'blue':
                median.set_color('orange')
            else:
                median.set_color('blue')
    plt.xticks(ticks=ticks, labels=tick_labels)
    plt.xlabel(plot_labels[0])
    plt.ylabel(plot_labels[1])
    plt.title(title)
    plt.legend([bplot1["boxes"][0], bplot2["boxes"][0]], labels, bbox_to_anchor=(1.14, 1.0), loc='upper right')
    plt.tight_layout()
    plt.savefig(out)
    plt.close()


def plot_vj_usage(xticks, sorted_values, dim, title, labels, outfile, tick_labels):
    plt.figure(figsize=dim)
    plt.xticks(xticks, labels=tick_labels, rotation=90)
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    xticks = numpy.arange(0.75, len(tick_labels) * 2, 2)
    means = []
    confidence_intervals = []
    for i in range(len(sorted_values[0])):
        mean, conf_int = calculate_confidence_intervals(sorted_values[0][i])
        means.append(mean)
        confidence_intervals.append(conf_int)
    plot_confidence_interval(xticks, sorted_values[0], means, confidence_intervals, 'blue', 'TdTKO')

    means = []
    confidence_intervals = []
    xticks = numpy.arange(1.25, (len(tick_labels)) * 2, 2)
    for i in range(len(sorted_values[1])):
        mean, conf_int = calculate_confidence_intervals(sorted_values[1][i])
        means.append(mean)
        confidence_intervals.append(conf_int)
    plot_confidence_interval(xticks, sorted_values[1], means, confidence_intervals, 'orange', 'WT')
    plt.subplots_adjust(bottom=0.15)
    plt.legend()
    plt.savefig(outfile)
    return


def plot_high_incidence_vj_usage(xticks, y_values, perc_no_d, dim, title, labels, outfile, tick_labels):
    plt.figure(figsize=dim)
    plt.xticks(xticks, labels=tick_labels, rotation=90)
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])

    plt.bar(xticks[0] - 0.2, y_values[0][0], color='blue', label='TdTKO', width=0.4)
    plt.bar(xticks[0] - 0.2, perc_no_d[0][0], color='black', width=0.4)
    plt.bar(xticks[0] + 0.2, y_values[1][0], color='orange', label='WT', width=0.4)
    plt.bar(xticks[0] + 0.2, perc_no_d[1][0], color='black', width=0.4)
    for i in range(1, len(xticks) - 1):
        plt.bar(xticks[i] - 0.2, y_values[0][i], color='blue', width=0.4)
        plt.bar(xticks[i] - 0.2, perc_no_d[0][i], color='black', width=0.4)
        plt.bar(xticks[i] + 0.2, y_values[1][i], color='orange', width=0.4)
        plt.bar(xticks[i] + 0.2, perc_no_d[1][i], color='black', width=0.4)

    plt.subplots_adjust(bottom=0.15)
    plt.legend()
    plt.savefig(outfile)
    return


def plot_vj_usage_with_without_d(xticks, y_values, dim, title, labels, outfile, tick_labels):
    plt.figure(figsize=dim)
    plt.xticks(xticks, labels=tick_labels, rotation=90)
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])

    plt.bar(xticks - 0.2, y_values[0], color='green', label='Inferred D segment = 0', width=0.4)
    plt.bar(xticks + 0.2, y_values[1], color='cyan', label='Inferred D segment > 0', width=0.4)

    plt.subplots_adjust(bottom=0.15)
    plt.legend()
    plt.savefig(outfile)
    return


def number_of_seq_per_incidence(x, y):

    plt.figure()
    plt.plot(x[0], y[0], label='TdTKO', color='blue')
    plt.scatter(x[0], y[0], color='blue', s=10)
    plt.plot(x[1], y[1], label='WT', color='orange')
    plt.scatter(x[1], y[1], color='orange', s=10)
    plt.yscale('log')
    plt.xticks(numpy.arange(0, 1.1, 0.1))
    plt.xlabel('Fraction of incidence')
    plt.ylabel('Percentage of sequences (%)')
    plt.title('Percentage of sequences per incidence')
    plt.legend()
    plt.savefig('C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\number_of_seq_per_incidence.png')


def plot_confidence_interval(xticks, values, means, confidence_intervals, colour, label, line_width=0.25):
    # mean = statistics.mean(values)
    # stdev = statistics.stdev(values)
    # confidence_interval = z * stdev / math.sqrt(len(values))
    # x = [[i + 1] * len(l[i]) for i in range(len(l))]
    for i in range(len(confidence_intervals)):
        left = xticks[i] - line_width / 2
        top = means[i] - confidence_intervals[i]
        right = xticks[i] + line_width / 2
        bottom = means[i] + confidence_intervals[i]
        plt.plot([xticks[i], xticks[i]], [top, bottom], color='black')
        plt.plot([left, right], [top, top], color='black')
        plt.plot([left, right], [bottom, bottom], color='black')
        plt.plot([left, right], [means[i], means[i]], color='green')
    xs = []
    ys = []
    for i in range(len(xticks)):
        try:
            x = [xticks[i]] * len(values[i])
            xs.extend(x)
            ys.extend(values[i])
        except IndexError as e:
            pass

    plt.scatter(xs, ys, color=colour, s=10, label=label)


def plot_deletions_conf_int(xticks, sorted_values,  title, labels, outfile, tick_labels, dim=(6.4, 4.8)):
    plt.figure(figsize=dim)
    # xticks = numpy.arange(1, len(tick_labels) + 1)
    plt.xticks(xticks, labels=tick_labels)
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    xticks = numpy.arange(0.75, len(tick_labels) * 2, 2)
    means = []
    confidence_intervals = []
    for i in range(len(sorted_values[0])):
        mean, conf_int = calculate_confidence_intervals(sorted_values[0][i])
        means.append(mean)
        confidence_intervals.append(conf_int)
    plot_confidence_interval(xticks, sorted_values[0], means, confidence_intervals, 'blue', 'TdTKO')

    means = []
    confidence_intervals = []
    xticks = numpy.arange(1.25, (len(tick_labels)) * 2, 2)
    for i in range(len(sorted_values[1])):
        mean, conf_int = calculate_confidence_intervals(sorted_values[1][i])
        means.append(mean)
        confidence_intervals.append(conf_int)
    plot_confidence_interval(xticks, sorted_values[1], means, confidence_intervals, 'orange', 'WT')
    plt.legend()
    plt.savefig(outfile)
    return


def plot_d_length_with_confidence_intervals(x_values:list, y_values, dim, title, labels, outfile):
    plt.figure(figsize=dim)
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    
    means = []
    confidence_intervals = []
    for i in range(len(x_values[0])):
        x = x_values[0][i]
        y = y_values[0][i]
        mean, conf_int = calculate_confidence_intervals(y)
        means.append(mean)
        confidence_intervals.append(conf_int)
    plot_confidence_interval(x, y, means, confidence_intervals, 'blue', 'TdTKO')
    plt.legend()
    plt.savefig(outfile)
    return


def plot_vj_usage_per_incidence(x_values, y_values, title, labels, outfile):
    plt.figure()
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.xticks(numpy.arange(0, 1.1, 0.1))
    plt.plot([n / len(x_values[0]) for n in x_values[0]], y_values[0], color='blue', label='TdTKO')
    plt.scatter([n / len(x_values[0]) for n in x_values[0]], y_values[0], color='blue', s=10)

    plt.plot([n / len(x_values[1]) for n in x_values[1]], y_values[1], color='orange', label='WT')
    plt.scatter([n / len(x_values[1]) for n in x_values[1]], y_values[1], color='orange', s=10)
    plt.legend()
    plt.savefig(outfile)
    plt.close()
    return


def plot_high_low_vj_usage(xticks, y_values, dim, title, labels, outfile, tick_labels):
    plt.figure(figsize=dim)
    plt.xticks(xticks, labels=tick_labels, rotation=90)
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.bar(xticks - 0.2, y_values[0], color='purple', label='High', edgecolor='black', width=0.4)
    plt.bar(xticks + 0.2, y_values[1], color='yellow', label='Low', edgecolor='black', width=0.4)
    plt.legend()
    plt.savefig(outfile)
    plt.close()
    return


def plot_pgens(data, pos, title, labels, outfile):
    plt.figure()
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.yscale('log')
    width = [0.05]*len(data)
    plt.boxplot(data, positions=pos, widths=width, showfliers=False, showmeans=True, meanline=True,
                vert=True, manage_ticks=False)
    plt.axhline(0, 0.5, 0.5, linestyle='dashed', color='green', label='Mean')
    plt.axhline(0, 0.5, 0.5, color='orange', label='Median')
    # plt.axhline()
    plt.legend()
    plt.savefig(outfile)
    plt.close()
    return


def plot_distribute_pgens(x_values, pos, title, labels, outfile):
    plt.figure()
    plt.title(title)
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    plt.yscale('log')
    plt.boxplot(x_values, positions=pos, showfliers=False, showmeans=True, meanline=True, vert=True)
    plt.axvline(0, 2.5e-07, 2.5e-07, linestyle='dashed', color='green', label='Mean')
    plt.axvline(0, 2.5e-07, 2.5e-07, color='orange', label='Median')
    plt.axhline(1.0e-40, 0, 0)
    plt.legend()
    plt.savefig(outfile)
    plt.close()
    return
