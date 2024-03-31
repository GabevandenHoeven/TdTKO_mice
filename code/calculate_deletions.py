import csv
from plots import plot_deletions, plot_d_deletions
import numpy


def get_vj_deletions(filename):
    """

    :param filename:
    :return:
    """
    v_del = {
        # deletion length: occurrences
    }
    j_del = {
        # deletion length: occurrences
    }
    total = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)
        for line in reader:
            total += 1
            try:
                v_del[int(line[header.index('V.length.deleted')])] += 1
            except KeyError:
                v_del.update({int(line[header.index('V.length.deleted')]): 1})
            try:
                j_del[int(line[header.index('J.length.deleted')])] += 1
            except KeyError:
                j_del.update({int(line[header.index('J.length.deleted')]): 1})
    return v_del, j_del, total


def calculate_d_deletions(filename):
    """

    :param filename:
    :return:
    """
    d_del = {}
    unidentified = 0
    d1_len = 12
    d2_len = 14
    total = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)
        for line in reader:
            total += 1
            d = line[header.index('D.used')]
            if 'T' in d or 'GGGGGG' in d:
                try:
                    deleted = d2_len - int(line[header.index('D.length.used')])
                    d_del[deleted] += 1
                except KeyError:
                    d_del.update({deleted: 1})
            elif 'CA' in d:
                try:
                    deleted = d1_len - int(line[header.index('D.length.used')])
                    d_del[deleted] += 1
                except KeyError:
                    d_del.update({deleted: 1})
            else:
                unidentified += 1
    return d_del, unidentified, total


if __name__ == '__main__':
    x = []
    y = []
    file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv'
    vdel_tdt, jdel_tdt, total_tdt = get_vj_deletions(file)
    x.append(sorted(vdel_tdt.keys()))
    y.append([vdel_tdt[e] / total_tdt * 100 for e in sorted(vdel_tdt.keys())])

    file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
    vdel_nor, jdel_nor, total_nor = get_vj_deletions(file)
    x.append(sorted(vdel_nor.keys()))
    y.append([vdel_nor[e] / total_nor * 100 for e in sorted(vdel_nor.keys())])

    out = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\deletion_size_V_TdTKO-Normal.png'
    title = 'Percentage of sequences containing deletions of V gene \nfor TdTKO and Normal'
    ticks = list(numpy.arange(14))
    plot_deletions(x, y, ['TdTKO', 'Normal'], ticks, title, out)

    x = []
    y = []

    x.append(sorted(jdel_tdt.keys()))
    y.append([jdel_tdt[e] / total_tdt * 100 for e in sorted(jdel_tdt.keys())])

    x.append(sorted(jdel_nor.keys()))
    y.append([jdel_nor[e] / total_nor * 100 for e in sorted(jdel_nor.keys())])
    out = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\deletion_size_J_TdTKO-Normal.png'
    title = 'Percentage of sequences containing deletions of J gene \nfor TdTKO and Normal'
    ticks = list(numpy.arange(14))
    plot_deletions(x, y, ['TdTKO', 'Normal'], ticks, title, out)

    # ----
    x = []
    y = []
    file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv'
    d_del_tdt, unidentified_tdt, total_tdt = calculate_d_deletions(file)
    x.append(sorted(d_del_tdt.keys()))
    x[-1].append(x[-1][-1] + 2)
    y.append([d_del_tdt[e] / total_tdt * 100 for e in sorted(d_del_tdt.keys())])
    y[-1].append(unidentified_tdt / total_tdt * 100)

    file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
    d_del_nor, unidentified_nor, total_nor = calculate_d_deletions(file)
    x.append(sorted(d_del_nor.keys()))
    x[-1].append(x[-1][-1] + 2)
    y.append([d_del_nor[e] / total_nor * 100 for e in sorted(d_del_tdt.keys())])
    y[-1].append(unidentified_nor / total_nor * 100)

    out = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\deletion_size_D_TdTKO-Normal.png'
    title = 'Percentage of sequences containing deletions of D gene \nfor TdTKO and Normal'
    plot_d_deletions(x, y, ['TdTKO', 'Normal'], title, out)
