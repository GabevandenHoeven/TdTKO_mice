import csv
from code_main import get_vdj_lengths, read_rtcr_refs


def filter_tcrb_data(rtcr_ref: dict, tcrb_data_filenames, threshold, max_sequence_size, flt_strain, delim, out_filename):
    """Reads the datafile with the nucleotide sequences and used vdj genes
    :param rtcr_ref: dict - reference sequences of TCR genes.
    :param tcrb_data_filenames: list - List with paths to the files containing the TCR data.
    :param threshold: int - Threshold for a minimum amount of supporting reads. More reads than this threshold are
    required to be included in the filtered data.
    :param max_sequence_size: int - The limit junction sequence length from which sequences are not included in
    the filtered data.
    :param flt_strain: str - What strain you want in the filtered data. (TdT-/-, C57BL/6Â )
    :param delim: str - The separation symbol for the data file.
    :param out_filename: str - Path to a new file in which the filtered data is stored.
    :return:
    """
    outlines = [
        "Column1\tMouse_ID\tNo_reads\tV_used\tD_used\tD_length\tJ_used\tV_deleted\tJ_deleted\t"
        "insertion_length\t""V_gene\tJ_gene\tV_end\tJ_start\tmin_phred\tcdr3\tphenotype\tstrain\n"]
    column = 0
    for file in tcrb_data_filenames:
        with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=delim)
            header = next(reader)
            for row in reader:
                column += 1
                reads, v_gene, j_gene, cdr3nt_seq, strain = row[header.index('Number.of.reads')], \
                    row[header.index('V.gene')], row[header.index('J.gene')], \
                    row[header.index('Junction.nucleotide.sequence')], row[header.index('strain')]
                line = get_vdj_lengths([v_gene, j_gene, cdr3nt_seq], rtcr_ref)
                if int(reads) > threshold and len(cdr3nt_seq) < max_sequence_size and strain == flt_strain:
                    v_end, j_start, phred, phenotype = row[header.index('V.gene.end.position')], \
                        row[header.index('J.gene.start.position')], row[header.index('Minimum.Phred')], \
                        row[header.index('phenotype')]
                    mouse_id = row[header.index('mouse')]

                    out_line = [str(column), mouse_id, reads]
                    out_line.extend(line)
                    out_line.extend([v_gene, j_gene, v_end, j_start, phred, cdr3nt_seq, phenotype, strain])
                    out_line = "\t".join(e for e in out_line) + "\n"
                    outlines.append(out_line)
    with open(out_filename, 'w') as outfile:
        outfile.writelines(outlines)
    return


if __name__ == '__main__':
    files = [
        # "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\tdt.csv"

        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\GSM6893351_Mandl-01-122016.tsv",
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\GSM6893356_Mandl-20-042016.tsv",
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\GSM6893359_Mandl-25-112016.tsv",
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\GSM6893360_Mandl-4-052016.tsv",
        "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\GSM6893362_Mandl-AF-human-07082016.tsv"

        # "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893369_RP-Mandl-28-M65&M66.tsv",
        # "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893370_RP-Mandl-30-M67&M68.tsv",
        # "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893365_RP-Mandl-05-M69&M70.tsv",
        # "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893366_RP-Mandl-06-M71&M72.tsv",
        # "C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\GSM6893367_RP-Mandl-07-M73&M74.tsv"
    ]
    refs = read_rtcr_refs()
    # new_filename = f'..\\data_files\\B6\\filtered_data\\filtered_data_Normal.tsv'
    new_filename = f'..\\data_files\\TdTKo\\filtered_data\\filtered_data_TdTKO.tsv'
    b6 = 'C57BL/6Â '
    tdt = 'TdT-/-'
    filter_tcrb_data(refs, files, 1, 64, tdt, '\t', new_filename)
