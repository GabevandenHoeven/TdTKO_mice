import csv

from code.code_main import get_vdj_lengths


def get_fractions(filename, reference_file, delim, threshold):
    """

    :param filename:
    :param reference_file:
    :param delim:
    :param threshold:
    :return:
    """
    total_above_threshold = 0
    insert_above_threshold = 0
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader)
        for row in reader:
            if row[header.index('Number.of.reads')] > threshold:
                total_above_threshold += 1
                v_, j_, crd3_ = row[header.index('V.gene')], row[header.index('J.gene')], \
                    row[header.index('Junction.nucleotide.sequence')]
                insertion_length = get_vdj_lengths([v_, j_, crd3_], reference_file)[6]
                if insertion_length > 0:
                    insert_above_threshold += 1
    return total_above_threshold, insert_above_threshold
