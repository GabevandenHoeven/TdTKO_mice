import csv
from utils import get_vdj_lengths


def read_tcrb_data(rtcr_ref: dict, tcrb_data_filename, data_suffix, threshold, delim):
    """Reads the datafile with the nucleotide sequences and used vdj genes
    :param rtcr_ref: dict - reference sequences of TCR genes.
    :param tcrb_data_filename: str - Path to the data containing the TCR data.
    :param data_suffix: str - suffix of the data data used to write results to a new data.
    :param threshold: int - Threshold for a minimum amount of supporting reads.
    :return:
    """

    with open(tcrb_data_filename, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=delim)
        header = next(reader)
        outfile = [
            "Column1\tMouse_ID\tNo_reads\tV_used\tD_used\tD_length\tJ_used\tV_deleted\tJ_deleted\tinsertion_length\tV_gene\tJ_gene\tV_end\tJ_start\tmin_phred\tcdr3\tphenotype\tstrain\tseq_id\n"]
        outfile_ins = [
            "Column1\tMouse_ID\tNo_reads\tV_used\tD_used\tD_length\tJ_used\tV_deleted\tJ_deleted\tinsertion_length\tV_gene\tJ_gene\tV_end\tJ_start\tmin_phred\tcdr3\tphenotype\tstrain\tseq_id\n"]
        outfile_below_threshold = [
            "Column1\tMouse_ID\tNo_reads\tV_used\tD_used\tD_length\tJ_used\tV_deleted\tJ_deleted\tinsertion_length\tV_gene\tJ_gene\tV_end\tJ_start\tmin_phred\tcdr3\tphenotype\tstrain\tseq_id\n"]
        supporting_read_count = {
            # number of reads: count of sequences
        }
        for row in reader:
            column, reads, v_gene, j_gene, cdr3nt_seq = row[0], \
                row[header.index('Number.of.reads')], row[header.index('V.gene')], row[header.index('J.gene')], \
                row[header.index('Junction.nucleotide.sequence')]
            try:
                mouse_id = row[header.index('mouse')]
            except ValueError:
                mouse_id = row[header.index('dirname')]

            line = get_vdj_lengths([v_gene, j_gene, cdr3nt_seq], rtcr_ref)
            out_line = [str(column), mouse_id, reads]
            out_line.extend(line)
            v_end, j_start, phred, phenotype, strain = row[header.index('V.gene.end.position')], \
                row[header.index('J.gene.start.position')], row[header.index('Minimum.Phred')], \
                row[header.index('phenotype')], row[header.index('strain')]
            out_line.extend([v_gene, j_gene, v_end, j_start, phred, cdr3nt_seq, phenotype, strain])
            out_line = "\t".join(e for e in out_line) + "\n"
            try:
                out_line = out_line.replace('\n', f'\t{row[header.index("Nucleotide.seq.id")]}\n')
            except ValueError:
                pass
            if int(reads) < threshold:
                outfile_below_threshold.append(out_line)
            elif line[-1] == "0":
                outfile.append(out_line)
            else:
                outfile_ins.append(out_line)
            try:
                supporting_read_count[int(reads)] += 1
            except KeyError:
                supporting_read_count.update({int(reads): 1})

    # with open('count_reads.txt', 'w') as data:
    #     data.write('Number_of_reads: number_of_sequences\n')
    #     for i in sorted(supporting_read_count.keys()):
    #         data.write(f'{i}: {supporting_read_count[i]}\n')

    # # new_filename = f'data_files\\B6\\no_insertions_above_threshold\\{data_suffix}_results_above_{threshold}.tsv'
    # new_filename = f'data_files\\TdTKo\\no_insertions_above_threshold\\{data_suffix}_results_above_{threshold}.tsv'
    # with open(new_filename, 'w') as data:
    #     data.writelines(outfile)
    # # new_insertion_filename = f'data_files\\B6\\insertions_above_threshold\\{data_suffix}_insertions_above_{threshold}.tsv'
    # new_insertion_filename = f'data_files\\TdTKo\\insertions_above_threshold\\{data_suffix}_insertions_above_{threshold}.tsv'
    # with open(new_insertion_filename, 'w') as data:
    #     data.writelines(outfile_ins)
    # # below_threshold_filename = f'data_files\\B6\\below_threshold\\{data_suffix}_below_{threshold}.tsv'
    # below_threshold_filename = f'data_files\\TdTKo\\below_threshold\\{data_suffix}_below_{threshold}.tsv'
    # with open(below_threshold_filename, 'w') as data:
    #     data.writelines(outfile_below_threshold)

    return


def fraction_insertions(all_sequences_filenames: list, insertions_filenames: list):
    """Calculates the fraction of sequences found containing insertions having above a given threshold of supporting
    reads.
    :param all_sequences_filenames: list - The path to the data containing all sequences of the sample
    :param insertions_filenames: list -
    """
    total_lines_ins = 0
    total_lines_all = 0
    for all_sequences_filename in all_sequences_filenames:
        with open(all_sequences_filename, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)

    for insertions_filename in insertions_filenames:
        with open(insertions_filename, 'r') as f:
            lines = len(f.readlines()) - 1
            total_lines_ins += lines
            print(lines)
    fraction = total_lines_ins * 100 / total_lines_all
    print(f'fraction = {total_lines_ins} * 100 / {total_lines_all} = {fraction}')


def get_junction_length(fn, lengths: list, reads: list, seq_counts: dict, delim: str):
    """

    :return:
    """
    with open(fn, 'r') as file:
        reader = csv.reader(file, delimiter=delim)
        header = next(reader)
        for row in reader:
            lengths.append(len(row[header.index('Junction.nucleotide.sequence')]))
            reads.append(int(row[header.index('Number.of.reads')]))
            try:
                seq_counts[len(row[header.index('Junction.nucleotide.sequence')])] += 1
            except KeyError:
                seq_counts.update({len(row[header.index('Junction.nucleotide.sequence')]): 1})
    return lengths, reads, seq_counts


if __name__ == "__main__":
    x_points = []
    y_points = []
    labels = []
    averages = []
    file_list = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO.tsv',
        '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal.tsv',
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal.tsv'
    ]

    # for filename in file_list:
    #     sequence_counts = {}
    #     sequence_lengths = []
    #     read_counts = []
    #     sequence_lengths, read_counts, sequence_counts = get_junction_length(filename, sequence_lengths, read_counts,
    #                                                                          sequence_counts, delim='\t')
    #     labels.append(filename.split('_')[-1].rstrip('.tsv'))
    #     x = sorted(sequence_counts.keys())
    #     x_points.append(x)
    #     y = [sequence_counts[point] for point in x]
    #     y = [point * 100 / sum(y) for point in y]
    #     y_points.append(y)
    #     average = sum(sequence_lengths) / len(sequence_lengths)
    #     averages.append(average)
    #     print(f'The average is: {average}')

    # plot_sequence_length_read_count(sequence_lengths, read_counts, label='Normal')
    # plot_dist_junction_sequence_length(x_points, y_points, labels, averages)

    # threshold_list = [0, 5, 10, 15, 20, 25, 50]
    # refs = read_rtcr_refs()
    # for threshold in threshold_list:
    #     read_tcrb_data(refs, filename, file_suffix, threshold)
    #     fraction_insertions(filename, file_suffix, threshold)
    #     fraction_insertions(file_list, files1)

    # reconfigure_header_tdt_normal_data(filename)

