import csv
from plots import plot_deletions, plot_d_deletions, plot_deletions_conf_int
import numpy
from utils import get_unique_sequences_per_mouse_from_file


def get_vj_deletions(data):
    """

    :param data:
    :return:
    """
    mice = {}

    total = 0
    header = data[0]
    for line in data[1:]:
        total += 1
        mouse = line[header.index('Mouse')]
        v_del, j_del = int(line[header.index('V.length.deleted')]), int(line[header.index('J.length.deleted')])
        if mouse not in mice.keys():
            mice.update({mouse: [{}, {}, 0]})
        mice[mouse][2] += 1
        if v_del not in mice[mouse][0].keys():
            mice[mouse][0].update({v_del: 1})
        else:
            mice[mouse][0][v_del] += 1
        if j_del not in mice[mouse][1].keys():
            mice[mouse][1].update({j_del: 1})
        else:
            mice[mouse][1][j_del] += 1
    return mice, total


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
    files = [
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]
    v_list = []
    j_list = []
    for file in files:
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        mice, total = get_vj_deletions(filtered_data)
        vs = []
        js = []
        for m in mice.keys():
            seqs = mice[m][2]
            vs.append([(mice[m][0][v] / seqs * 100) for v in sorted(mice[m][0].keys())])
            js.append([(mice[m][1][j] / seqs * 100) for j in sorted(mice[m][1].keys())])
        v_list.append(vs)
        j_list.append(js)

    sorted_values = []
    for i in range(len(v_list)):
        sorted_values_i = []
        for ii in range(max([len(e) for e in v_list[i]])):
            x = [v_list[i][m][ii] for m in range(len(v_list[i]))]
            sorted_values_i.append(x)
        sorted_values.append(sorted_values_i)

    xticks = numpy.arange(1, (max(len(e) for e in sorted_values)) * 2, 2)
    x_labels = numpy.arange((max(len(e) for e in sorted_values)))
    plot_deletions_conf_int(xticks, sorted_values, (8, 10), 'Deletions in V segments with confidence intervals',
                            ('Deletion size (nt)', 'Percentage sequences (%)'),
                            '..\\img\\unique_seq_img\\V_del_confint_gen.png', x_labels)
    # plot_deletions_conf_int(xticks, sorted_values, (8, 10), 'V deletions confint exp',
    #                         ('Deletion size (nt)', 'Percentage sequences (%)'),
    #                         '..\\img\\V_del_confint.png', x_labels)

    sorted_values = []
    for i in range(len(j_list)):
        sorted_values_i = []
        for ii in range(max([len(e) for e in j_list[i]])):
            x = []
            for m in range(len(j_list[i])):
                if len(j_list[i][m]) > ii:
                    x.append(j_list[i][m][ii])
            # x2 = [j_list[i][m][ii] for m in range(len(j_list[i]))]
            sorted_values_i.append(x)
        sorted_values.append(sorted_values_i)

    xticks = numpy.arange(1, (max(len(e) for e in sorted_values)) * 2, 2)
    x_labels = numpy.arange((max(len(e) for e in sorted_values)))
    plot_deletions_conf_int(xticks, sorted_values, (8, 10), 'Deletions in J segments with confidence intervals',
                            ('Deletion size (nt)', 'Percentage sequences (%)'),
                            '..\\img\\unique_seq_img\\J_del_confint_gen.png', x_labels)
    # plot_deletions_conf_int(xticks, sorted_values, (8, 10), 'J deletions exp',
    #                         ('Deletion size (nt)', 'Percentage sequences (%)'),
    #                         '..\\img\\J_del_confint.png', x_labels)

    # ----------------------------------------------------------------------------------

    # file = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
    # # file = '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal.tsv'
    # vdel_nor, jdel_nor, total_nor = get_vj_deletions(file)
    # x.append(sorted(vdel_nor.keys()))
    # y.append([vdel_nor[e] / total_nor * 100 for e in sorted(vdel_nor.keys())])
    #
    # out = f'..\\img\\deletion_size_V_TdTKO-Normal2.png'
    # # out = f'..\\img\\Generated_deletion_size_V_TdTKO-Normal.png'
    # title = 'Percentage of sequences containing deletions of V gene \nfor TdTKO and Normal'
    # ticks = list(numpy.arange(14))
    # # plot_deletions(x, y, ['TdTKO', 'Normal'], ticks, title, out)
    #
    # x = []
    # y = []
    #
    # x.append(sorted(jdel_tdt.keys()))
    # y.append([jdel_tdt[e] / total_tdt * 100 for e in sorted(jdel_tdt.keys())])
    #
    # x.append(sorted(jdel_nor.keys()))
    # y.append([jdel_nor[e] / total_nor * 100 for e in sorted(jdel_nor.keys())])
    # out = f'..\\img\\deletion_size_J_TdTKO-Normal2.png'
    # # out = f'..\\img\\Generated_deletion_size_J_TdTKO-Normal.png'
    # title = 'Percentage of sequences containing deletions of J gene \nfor TdTKO and Normal'
    # ticks = list(numpy.arange(14))
    # # plot_deletions(x, y, ['TdTKO', 'Normal'], ticks, title, out)

    # ----
    # x = []
    # y = []
    # file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv'
    # file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\generated_filtered_data_TdTKO.tsv'
    # d_del_tdt, unidentified_tdt, total_tdt = calculate_d_deletions(file)
    # x.append(sorted(d_del_tdt.keys()))
    # x[-1].append(x[-1][-1] + 2)
    # y.append([d_del_tdt[e] / total_tdt * 100 for e in sorted(d_del_tdt.keys())])
    # y[-1].append(unidentified_tdt / total_tdt * 100)

    # file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
    # file = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\generated_filtered_data_Normal.tsv'
    # d_del_nor, unidentified_nor, total_nor = calculate_d_deletions(file)
    # x.append(sorted(d_del_nor.keys()))
    # x[-1].append(x[-1][-1] + 2)
    # y.append([d_del_nor[e] / total_nor * 100 for e in sorted(d_del_tdt.keys())])
    # y[-1].append(unidentified_nor / total_nor * 100)

    # out = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\deletion_size_D_TdTKO-Normal.png'
    # out = f'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\img\\Generated_deletion_size_D_TdTKO-Normal.png'
    # title = 'Percentage of sequences containing deletions of D gene \nfor TdTKO and Normal'
    # plot_d_deletions(x, y, ['TdTKO', 'Normal'], title, out)
