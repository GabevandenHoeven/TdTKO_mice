import numpy

from utils import get_unique_sequences_from_file


def get_insertions_per_d(data):
    ins_per_d = {}
    header = data[0]
    for line in data[1:]:
        insertions = int(line[header.index('Left.insertion.length')]) + \
                     int(line[header.index('Right.insertion.length')])
        d_length = int(line[header.index('D.length.used')])
        try:
            ins_per_d[d_length].append(insertions)
        except KeyError:
            ins_per_d.update({d_length: [insertions]})
    return ins_per_d


if __name__ == '__main__':

    filename = '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v2.tsv'
    filtered_data = get_unique_sequences_from_file(filename)
    insertion_data = get_insertions_per_d(filtered_data)

    print('Mean insertion lengths per inferred D segment length')
    for d in sorted(insertion_data.keys()):
        print(f'{d}: {sum(insertion_data[d]) / len(insertion_data[d])}',
              len(insertion_data[d]), numpy.median(insertion_data[d]))
