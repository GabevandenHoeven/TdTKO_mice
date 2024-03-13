import csv
from plots import plot_deletions


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
    title = 'Percentage of sequences containing deletions of V gene for TdTKO and Normal'
    plot_deletions(x, y, ['TdTKO', 'Normal'], title, out)

    x = []
    y = []

    x.append(sorted(jdel_tdt.keys()))
    y.append([jdel_tdt[e] / total_tdt * 100 for e in sorted(jdel_tdt.keys())])

    x.append(sorted(jdel_nor.keys()))
    y.append([jdel_nor[e] / total_nor * 100 for e in sorted(jdel_nor.keys())])
    out = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\deletion_size_J_TdTKO-Normal.png'
    title = 'Percentage of sequences containing deletions of J gene for TdTKO and Normal'
    plot_deletions(x, y, ['TdTKO', 'Normal'], title, out)
