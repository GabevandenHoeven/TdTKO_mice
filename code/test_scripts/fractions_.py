import csv

# from  import get_vdj_lengths, read_rtcr_refs


def get_fractions(filename, reference_file, delim, threshold, total_above_threshold, insert_above_threshold):
    """

    :param filename:
    :param reference_file:
    :param delim:
    :param threshold:
    :param total_above_threshold:
    :param insert_above_threshold:
    :return:
    """
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader)
        for row in reader:
            if int(row[header.index('Number.of.reads')]) > threshold:
                total_above_threshold += 1
                v_, j_, crd3_ = row[header.index('V.gene')], row[header.index('J.gene')], \
                    row[header.index('Junction.nucleotide.sequence')]
                insertion_length = get_vdj_lengths([v_, j_, crd3_], reference_file)[6]
                if int(insertion_length) > 0:
                    insert_above_threshold += 1
    return total_above_threshold, insert_above_threshold


if __name__ == '__main__':
    refs = read_rtcr_refs()
    threshold_list = [0, 1, 2, 5, 10, 15, 20, 25, 50]
    files = [
        # 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\tdt.csv'
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893369_RP-Mandl-28-M65&M66.tsv",
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893370_RP-Mandl-30-M67&M68.tsv",
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893365_RP-Mandl-05-M69&M70.tsv",
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893366_RP-Mandl-06-M71&M72.tsv",
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893367_RP-Mandl-07-M73&M74.tsv"
    ]
    # Use lists to get all thresholds in one go without having to read the same files multiple times
    for threshold in threshold_list:
        above_t = 0
        ins_above_t = 0
        for file in files:
            above_t, ins_above_t = get_fractions(file, refs, '\t', threshold, above_t, ins_above_t)
        fract = ins_above_t * 100 / above_t
        print(f'fraction for threshold: {threshold} = {ins_above_t} * 100 / {above_t} = {fract}')
