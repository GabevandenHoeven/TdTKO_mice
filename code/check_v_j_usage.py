import csv
import matplotlib.pyplot as plt
import numpy
from plots import plot_vj_usage, plot_confidence_interval


def check_vj_usage(filename: str, delim='\t'):
    """Check what V and J genes are used in each file
    :param filename:
    :param delim: str
    :return:
    """
    mice = {}
    count_lines = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            count_lines += 1
            mouse = line[header.index('Mouse')]
            v_gene, j_gene = line[header.index('V.gene')], line[header.index('J.gene')]
            if mouse not in mice.keys():
                mice.update({mouse: [{}, {}, 0]})
            mice[mouse][2] += 1
            if v_gene not in mice[mouse][0].keys():
                mice[mouse][0].update({v_gene: 1})
            else:
                mice[mouse][0][v_gene] += 1
            if j_gene not in mice[mouse][1].keys():
                mice[mouse][1].update({j_gene: 1})
            else:
                mice[mouse][1][j_gene] += 1
    return mice, count_lines


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv',
        # '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        # '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]
    v_list = []
    j_list = []
    for file in files:
        mice, count = check_vj_usage(file)
        fn = file.split('\\')[-1]
        vs = []
        js = []
        for m in mice.keys():
            seqs = mice[m][2]
            vs.append([(mice[m][0][v] / seqs * 100) for v in sorted(mice[m][0].keys())])
            js.append([(mice[m][1][j] / seqs * 100) for j in sorted(mice[m][1].keys())])
        v_list.append(vs)
        j_list.append(js)
    v_labels = ['TRBV1*01', 'TRBV12-1*01', 'TRBV12-2*01', 'TRBV13-1*01', 'TRBV13-2*01', 'TRBV13-3*01', 'TRBV14*01',
                'TRBV15*01', 'TRBV16*01', 'TRBV17*01', 'TRBV19*01', 'TRBV2*01', 'TRBV20*01', 'TRBV23*01', 'TRBV24*01',
                'TRBV26*01', 'TRBV29*01', 'TRBV3*01', 'TRBV30*01', 'TRBV31*01', 'TRBV4*01', 'TRBV5*01']
    j_labels = ['TRBJ1-1*01', 'TRBJ1-2*01', 'TRBJ1-3*01', 'TRBJ1-4*01', 'TRBJ1-5*01', 'TRBJ2-1*01', 'TRBJ2-2*01',
                'TRBJ2-3*01', 'TRBJ2-4*01', 'TRBJ2-5*01', 'TRBJ2-7*01']
    sorted_values = []
    for i in range(len(v_list)):
        sorted_values_i = []
        for segment in range(len(v_list[i][0])):
            x = [v_list[i][m][segment] for m in range(len(v_list[i]))]
            sorted_values_i.append(x)
        sorted_values.append(sorted_values_i)
    xticks = numpy.arange(1, (len(v_labels)) * 2, 2)
    plot_vj_usage(xticks, sorted_values, (8, 10), 'V usage exp', ('V segments', 'Usage (%)'),
                  '..\\img\\V_usage_exp.png', v_labels)

    sorted_values = []
    for i in range(len(j_list)):
        sorted_values_i = []
        for segment in range(len(j_list[i][0])):
            x = [j_list[i][m][segment] for m in range(len(j_list[i]))]
            sorted_values_i.append(x)
        sorted_values.append(sorted_values_i)
    xticks = numpy.arange(1, (len(j_labels)) * 2, 2)
    plot_vj_usage(xticks, sorted_values, (8, 10), 'J usage exp', ('J segments', 'Usage (%)'),
                  '..\\img\\J_usage_exp.png', j_labels)
