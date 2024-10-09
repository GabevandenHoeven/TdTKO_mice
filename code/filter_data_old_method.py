import csv
from utils import get_vdj_lengths_old_method, read_rtcr_refs


def filter_tcrb_data(rtcr_ref: dict, tcrb_data_filenames, threshold, max_sequence_size, flt_strain, flt_pheno, delim,
                     out_filename):
    """Reads the datafile with the nucleotide sequences and used vdj genes
    :param rtcr_ref: dict - reference sequences of TCR genes.
    :param tcrb_data_filenames: list - List with paths to the files containing the TCR data.
    :param threshold: int - Threshold for a minimum amount of supporting reads. More reads than this threshold are
    required to be included in the filtered data.
    :param max_sequence_size: int - The limit junction sequence length from which sequences are not included in
    the filtered data.
    :param flt_strain: str - What strain you want in the filtered data. (TdT-/-, C57BL/6Â )
    :param flt_pheno: str - What phenotype cell to want in the filtered data. (CD4+, CD5hi, CD5lo)
    :param delim: str - The separation symbol for the datafile.
    :param out_filename: str - Path to a new file in which the filtered data is stored.
    :return:
    """
    outlines = [
        "Column1\tMouse\tNumber.of.reads\tV.length.used\tD.used\tD.length.used\tJ.length.used\tV.length.deleted\t"
        "J.length.deleted\tInsertions\tV.J.distance\tV.gene\tJ.gene\tV.gene.end.position\t"
        "J.gene.start.position\tMin.Phred\tJunction.nucleotide.sequence\tphenotype\tstrain\n"]
    column = 0
    for file in tcrb_data_filenames:
        with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=delim)
            header = next(reader)
            for row in reader:
                reads, v_gene, j_gene, cdr3nt_seq, strain, phenotype, phred = row[header.index('Number.of.reads')], \
                    row[header.index('V.gene')], row[header.index('J.gene')], \
                    row[header.index('Junction.nucleotide.sequence')], row[header.index('strain')], \
                    row[header.index('phenotype')], row[header.index('Minimum.Phred')]
                line = get_vdj_lengths_old_method([v_gene, j_gene, cdr3nt_seq], rtcr_ref)
                v_end, j_start = line[-2], line[-1]
                line.pop(-1), line.pop(-1)
                if int(reads) > threshold and len(cdr3nt_seq) < max_sequence_size and strain == flt_strain and \
                        phenotype == flt_pheno:
                    column += 1
                    mouse_id = row[header.index('mouse')]
                    v_j_distance = str(int(line[2]) + int(line[6]))
                    out_line = [str(column), mouse_id, reads]
                    out_line.extend(line)
                    out_line.extend([v_j_distance, v_gene, j_gene, v_end, j_start, phred,
                                     cdr3nt_seq, phenotype, strain])
                    out_line = "\t".join(e for e in out_line) + "\n"
                    outlines.append(out_line)
    with open(out_filename, 'w') as outfile:
        outfile.writelines(outlines)
    return


def filter_tcrb_generated_data(rtcr_ref: dict, tcrb_data_filenames, max_sequence_size, flt_strain,
                               delim, out_filename):
    """Reads the datafile with the nucleotide sequences and used vdj genes
    :param rtcr_ref: dict - reference sequences of TCR genes.
    :param tcrb_data_filenames: list - List with paths to the files containing the TCR data.
    :param max_sequence_size: int - The limit junction sequence length from which sequences are not included in
    the filtered data.
    :param flt_strain: str - What strain you want in the filtered data. (TdT-/-, C57BL/6Â )
    :param delim: str - The separation symbol for the datafile.
    :param out_filename: str - Path to a new file in which the filtered data is stored.
    :return:
    """
    outlines = [
        "Column1\tMouse\tV.length.used\tD.used\tD.length.used\tJ.length.used\tV.length.deleted\t"
        "J.length.deleted\tInsertions\tV.J.distance\tD.used.data\tV.gene\tD.gene\tJ.gene\tV.gene.end.position\t"
        "J.gene.start.position\tJunction.nucleotide.sequence\tphenotype\tstrain\n"]
    column = 0
    for file in tcrb_data_filenames:
        with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=delim)
            header = next(reader)
            for row in reader:
                v_gene, j_gene, d_gene, cdr3nt_seq, strain, phenotype, d_used_data = row[header.index('V.gene')], \
                    row[header.index('J.gene')], row[header.index('D.gene')], \
                    row[header.index('Junction.nucleotide.sequence')], row[header.index('strain')], \
                    row[header.index('phenotype')], row[header.index('D.used')]
                line = get_vdj_lengths_old_method([v_gene, j_gene, cdr3nt_seq], rtcr_ref)
                if len(cdr3nt_seq) < max_sequence_size and strain == flt_strain:
                    column += 1
                    v_end, j_start = line[-2], line[-1]
                    line.pop(-1), line.pop(-1)
                    mouse = row[header.index('mouse')]
                    v_j_distance = str(int(line[6]) + int(line[2]) + int(line[7]))
                    out_line = [str(column), mouse]
                    out_line.extend(line)
                    out_line.extend([v_j_distance, d_used_data, v_gene, d_gene, j_gene, v_end, j_start,
                                     cdr3nt_seq, phenotype, strain])
                    out_line = "\t".join(out_line) + "\n"
                    outlines.append(out_line)
    with open(out_filename, 'w') as outfile:
        outfile.writelines(outlines)
        return


if __name__ == '__main__':
    refs = read_rtcr_refs()

    files = [
        "..\\data_files\\Normal\\GSM6893369_RP-Mandl-28-M65&M66.tsv",
        "..\\data_files\\Normal\\GSM6893370_RP-Mandl-30-M67&M68.tsv",
        "..\\data_files\\Normal\\GSM6893365_RP-Mandl-05-M69&M70.tsv",
        "..\\data_files\\Normal\\GSM6893366_RP-Mandl-06-M71&M72.tsv",
        "..\\data_files\\Normal\\GSM6893367_RP-Mandl-07-M73&M74.tsv"
    ]

    new_filename = f'..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal_v1_test.tsv'
    b6 = 'C57BL/6Â '
    pheno = 'CD4+'
    filter_tcrb_data(refs, files, 1, 64, b6, pheno, '\t', new_filename)

    files = [
        "..\\data_files\\TdTKO\\GSM6893351_Mandl-01-122016.tsv",
        "..\\data_files\\TdTKO\\GSM6893356_Mandl-20-042016.tsv",
        "..\\data_files\\TdTKO\\GSM6893359_Mandl-25-112016.tsv",
        "..\\data_files\\TdTKO\\GSM6893360_Mandl-4-052016.tsv",
        "..\\data_files\\TdTKO\\GSM6893362_Mandl-AF-human-07082016.tsv"
    ]
    new_filename = f'..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO_v1_test.tsv'
    tdt = 'TdT-/-'
    filter_tcrb_data(refs, files, 1, 64, tdt, pheno, '\t', new_filename)


