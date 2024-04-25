from olga import load_model, sequence_generation as seq_gen


def convert_olga_to_imgt_id(olga_id: int, gene: str):
    """Converts the OLGA identifiers of V, D, and J segments to the general identifiers used on IMGT.

    :param olga_id: int - number representing the gene segment.
    :param gene: str - 'V', 'D', or 'J' used to look up the ID in the correct dictionary.
    :return:
    """
    # I checked the V gene segments, these are the right indices for these gene segments. Some gene segments have the
    #  same sequences they are mentioned in the comments below but these are the right indices
    v_gene_segments = {
        0: 'TRBV1*01', 1: 'TRBV2*01', 2: 'TRBV3*01', 3: 'TRBV4*01', 4: 'TRBV5*01',
        # index TRBV2*01 and TRBV5*01 have the same sequence
        11: 'TRBV12-1*01', 12: 'TRBV13-1*01', 13: 'TRBV12-2*01', 14: 'TRBV13-2*01', 16: 'TRBV13-3*01', 17: 'TRBV14*01',
        # TRBV12-1*01 and TRBV12-2*01 have the same sequence
        # index TRBV13-1*01 and TRBV13-3*01 have the same sequence
        18: 'TRBV15*01', 19: 'TRBV16*01', 20: 'TRBV17*01', 22: 'TRBV19*01', 23: 'TRBV20*01', 26: 'TRBV23*01',
        27: 'TRBV24*01', 29: 'TRBV26*01', 32: 'TRBV29*01', 33: 'TRBV30*01', 34: 'TRBV31*01'
        # TRBV31*01 sequence in OLGA default_model mouse_T_beta doesn't match the IMGT sequence
    }
    j_gene_segments = {
        0: 'TRBJ1-1*01', 1: 'TRBJ1-2*01', 2: 'TRBJ1-3*01', 3: 'TRBJ1-4*01', 4: 'TRBJ1-5*01',
        7: 'TRBJ2-1*01', 8: 'TRBJ2-2*01', 9: 'TRBJ2-3*01', 10: 'TRBJ2-4*01', 11: 'TRBJ2-5*01', 13: 'TRBJ2-7*01'
    }
    d_gene_segments = {
        0: 'TRBD1*01', 1: 'TRBD2*01'
    }
    match gene:
        case 'V':
            return v_gene_segments[olga_id]
        case 'J':
            return j_gene_segments[olga_id]
        case 'D':
            return d_gene_segments[olga_id]


if __name__ == '__main__':
    params_filename = 'models\\Mus+musculus\\TRB\\models\\model_params.txt'
    marginals_filename = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\model_marginals.txt'
    v_anchor_pos_file = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\V_gene_CDR3_anchors.csv'
    j_anchor_pos_file = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\J_gene_CDR3_anchors.csv'
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_filename, v_anchor_pos_file, j_anchor_pos_file)
    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_filename)
    seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)

    n = 1000000
    filename = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\generated_data_TdTKO.tsv'
    with open(filename, 'w') as outfile:
        outfile.write('Column1\tV.gene\tD.gene\tJ.gene\tJunction.nucleotide.sequence\tphenotype\tstrain\n')
    with open(filename, 'a') as outfile:
        for i in range(1, n + 1):
            sequence, aa, v, j, d = seq_gen_model.gen_rnd_prod_noins_CDR3()
            v = convert_olga_to_imgt_id(v, 'V')
            d = convert_olga_to_imgt_id(d, 'D')
            j = convert_olga_to_imgt_id(j, 'J')
            outfile.write(f'{i}\t{v}\t{d}\t{j}\t{sequence}\tGenerated\tTdT-/-\n')

    filename = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\generated_data_Normal.tsv'
    with open(filename, 'w') as outfile:
        outfile.write('Column1\tV.gene\tD.gene\tJ.gene\tJunction.nucleotide.sequence\tphenotype\tstrain\n')
    with open(filename, 'a') as outfile:
        for i in range(n):
            sequence, aa, v, j, d = seq_gen_model.gen_rnd_prod_CDR3()
            v = convert_olga_to_imgt_id(v, 'V')
            d = convert_olga_to_imgt_id(d, 'D')
            j = convert_olga_to_imgt_id(j, 'J')
            outfile.write(f'{i}\t{v}\t{d}\t{j}\t{sequence}\tGenerated\tC57BL/6Â \n')
