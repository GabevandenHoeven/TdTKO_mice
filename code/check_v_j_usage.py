import numpy
from plots import plot_vj_usage
from get_olga_vj_usage import get_olga_vj_usage
from utils import get_unique_sequences_per_mouse_from_file


def check_vj_usage(data):
    """Check what V and J genes are used in a datafile. Returns a hashmap with V and J segment usage per mouse.
    :param data: nested list - A nested list of lines from a datafile. The first line is a header.
    :return:
    """
    mice_hashmap = {}
    count_lines = 0
    header = data[0]
    for line in data[1:]:
        count_lines += 1
        mouse = line[header.index('Mouse')]
        v_gene, j_gene = line[header.index('V.gene')], line[header.index('J.gene')]
        if mouse not in mice_hashmap.keys():
            mice_hashmap.update({mouse: [{}, {}, 0]})
        mice_hashmap[mouse][2] += 1
        if v_gene not in mice_hashmap[mouse][0].keys():
            mice_hashmap[mouse][0].update({v_gene: 1})
        else:
            mice_hashmap[mouse][0][v_gene] += 1
        if j_gene not in mice_hashmap[mouse][1].keys():
            mice_hashmap[mouse][1].update({j_gene: 1})
        else:
            mice_hashmap[mouse][1][j_gene] += 1
    return mice_hashmap, count_lines


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv',
    ]
    v_list = []
    j_list = []
    for file in files:
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        mice, count = check_vj_usage(filtered_data)
        fn = file.split('\\')[-1]
        vs = []
        js = []
        for m in mice.keys():
            seqs = mice[m][2]
            vs.append([(mice[m][0][v] / seqs * 100) for v in sorted(mice[m][0].keys())])
            js.append([(mice[m][1][j] / seqs * 100) for j in sorted(mice[m][1].keys())])
        v_list.append(vs)
        j_list.append(js)
    olga_v, olga_j = get_olga_vj_usage()
    v_labels = ['TRBV1', 'TRBV12-1', 'TRBV12-2', 'TRBV13-1', 'TRBV13-2', 'TRBV13-3', 'TRBV14',
                'TRBV15', 'TRBV16', 'TRBV17', 'TRBV19', 'TRBV2', 'TRBV20', 'TRBV23', 'TRBV24',
                'TRBV26', 'TRBV29', 'TRBV3', 'TRBV30', 'TRBV31', 'TRBV4', 'TRBV5']
    j_labels = ['TRBJ1-1', 'TRBJ1-2', 'TRBJ1-3', 'TRBJ1-4', 'TRBJ1-5', 'TRBJ2-1', 'TRBJ2-2',
                'TRBJ2-3', 'TRBJ2-4', 'TRBJ2-5', 'TRBJ2-7']
    sorted_values = []
    for i in range(len(v_list)):
        sorted_values_i = []
        for segment in range(len(v_list[i][0])):
            x = [v_list[i][m][segment] for m in range(len(v_list[i]))]
            sorted_values_i.append(x)
        sorted_values.append(sorted_values_i)
    xticks = numpy.arange(1, (len(v_labels)) * 2, 2)
    plot_vj_usage(xticks, sorted_values, (8, 10), 'V usage TdTKO vs WT', ('V segments', 'Usage (%)'),
                  '..\\img\\V_usage.png', v_labels)

    sorted_values = []
    for i in range(len(j_list)):
        sorted_values_i = []
        for segment in range(len(j_list[i][0])):
            x = [j_list[i][m][segment] for m in range(len(j_list[i]))]
            sorted_values_i.append(x)
        sorted_values.append(sorted_values_i)
    xticks = numpy.arange(1, (len(j_labels)) * 2, 2)
    plot_vj_usage(xticks, sorted_values, (8, 10), 'J usage TdTKO vs WT', ('J segments', 'Usage (%)'),
                  '..\\img\\J_usage.png', j_labels)
