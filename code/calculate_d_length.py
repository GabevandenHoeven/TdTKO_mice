from utils import get_unique_sequences_from_file


def get_d_length_and_fraction_with_d(data, threshold):
    """Calculates the mean D segment length and the fraction of sequences that have a D segment, for a given threshold
    of considering an inferred D segment true.

    :param data: nested list - A nested list of lines from a datafile. The first line is a header.
    :param threshold: int - A given length of a D segment as a threshold for considering it a true D segment.
    :return:
    """
    avg_d_length = []
    n_lines = 0
    n_d_seg = 0
    header = data[0]
    for line in data[1:]:
        n_lines += 1
        if int(line[header.index('D.length.used')]) >= threshold:
            avg_d_length.append(int(line[header.index('D.length.used')]))
            n_d_seg += 1
    fract_d = n_d_seg * 100 / n_lines
    print(f'The average D length is: {sum(avg_d_length) / len(avg_d_length)}')
    print(f'The fraction of sequences that have a D segment is: {n_d_seg} * 100 / {n_lines} = {fract_d}')


if __name__ == '__main__':
    files = [
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv',
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v2.tsv',
        # '..\\data_files\\WT\\gen_data_mess\\filtered_data_gen_WT_v2.tsv',
        # '..\\data_files\\TdTKO\\gen_data_mess\\filtered_data_gen_TdTKO_v2.tsv'
    ]
    for file in files:
        filtered_data = get_unique_sequences_from_file(file)
        get_d_length_and_fraction_with_d(filtered_data, 1)
