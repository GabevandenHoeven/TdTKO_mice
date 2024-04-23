from olga import load_model, sequence_generation as seq_gen


def convert_olga_to_imgt_id(olga_id: int, gene: str):
    """Converts the OLGA identifiers of V, D, and J segments to the general identifiers used on IMGT.

    :param olga_id: int - number representing the gene segment.
    :param gene: str - 'V', 'D', or 'J' used to look up the ID in the correct dictionary.
    :return:
    """
    v_gene_segments = {
        0: '', 1: '', 2: '', 3: '', 4: '',
        11: '', 12: '', 13: '', 14: '', 16: '', 17: '', 18: '', 19: '',
        20: '', 22: '', 23: '', 26: '', 27: '', 29: '',
        32: '', 33: '', 34: ''
    }
    j_gene_segments = {
        0: '', 1: '', 2: '', 3: '', 4: '', 7: '', 8: '', 9: '',
        10: '', 11: '', 13: ''
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
    params_filename = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\model_params.txt'
    marginals_filename = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\model_marginals.txt'
    v_anchor_pos_file = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\V_gene_CDR3_anchors.csv'
    j_anchor_pos_file = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\J_gene_CDR3_anchors.csv'
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_filename, v_anchor_pos_file, j_anchor_pos_file)
    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_filename)
    seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)

    n = 500
    vs = []
    ds = []
    js = []
    filename = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\TdTKo\\filtered_data\\generated_data_TdTKO.tsv'
    with open(filename, 'w') as outfile:
        outfile.write('Column1\tV.gene\tD.gene\tJ.gene\tJunction.nucleotide.sequence\tphenotype\tstrain\n')
    with open(filename, 'a') as outfile:
        for i in range(1, n + 1):
            sequence, aa, v, j, d = seq_gen_model.gen_rnd_prod_CDR3_noins()
            # v = convert_olga_to_imgt_id(v, 'V')
            # d = convert_olga_to_imgt_id(d, 'D')
            # j = convert_olga_to_imgt_id(j, 'J')
            # vs.append(v)
            # ds.append(d)
            # js.append(j)
            # outfile.write(f'{i}\t{v}\t{d}\t{j}\t{sequence}\tGenerated\tTdTKO\n')
            # print(f'seq: {sequence}\nv: {v}\nd: {d}\nj: {j}')
        # print(set(vs), len(set(vs)))
        # print(set(ds), len(set(ds)))
        # print(set(js), len(set(js)))

    # filename = 'C:\\Users\\gabev\\PycharmProjects\\MRP_TdTKO_mice\\data_files\\B6\\filtered_data\\generated_data_Normal.tsv'
    # for i in range(n):
    #     sequence, aa, v, j, d = seq_gen_model.gen_rnd_prod_CDR3()
