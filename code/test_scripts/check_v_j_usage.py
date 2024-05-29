import csv


def check_vj_usage(filename: str, delim='\t'):
    """Check what V and J genes are used in each file
    :param filename:
    :param delim: str
    :return:
    """
    v_usage = {}
    j_usage = {}
    count_lines = 0
    with open(filename, 'r') as infile:
        reader = csv.reader(infile, delimiter=delim)
        header = next(reader)
        for line in reader:
            count_lines += 1
            v_gene, j_gene = line[header.index('V.gene')], line[header.index('J.gene')]
            try:
                v_usage[v_gene] += 1
            except KeyError:
                v_usage.update({v_gene: 1})
            try:
                j_usage[j_gene] += 1
            except KeyError:
                j_usage.update({j_gene: 1})
    return v_usage, j_usage, count_lines


if __name__ == '__main__':
    files = [
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv',
        '..\\..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO_v2.tsv',
        '..\\..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal_v2.tsv'
    ]
    for file in files:
        v, j, count = check_vj_usage(file)
        fn = file.split('\\')[-1]
        for i in v.keys():
            v[i] = v[i] / count * 100
        for i in j.keys():
            j[i] = j[i] / count * 100
        print(f'VJ usage in file: {fn}\nV: {[(i, v[i]) for i in sorted(v.keys())]}\nJ: {[(i, j[i]) for i in sorted(j.keys())]}')
