import threading
import time

from utils import animate
from olga import load_model, generation_probability as pgen_mod
import csv


def load_olga_mouse_trb_pgen_model():
    params_filename = 'models\\Mus+musculus\\TRB\\models\\model_params.txt'
    marginals_filename = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\model_marginals.txt'
    v_anchor_pos_file = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\V_gene_CDR3_anchors.csv'
    j_anchor_pos_file = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\J_gene_CDR3_anchors.csv'
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_filename, v_anchor_pos_file, j_anchor_pos_file)
    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_filename)
    pgen_model = pgen_mod.GenerationProbabilityVDJ(generative_model, genomic_data)
    return pgen_model


def calculate_pgen(filename, pgen):
    with open(filename, 'r') as in_file:
        reader = csv.reader(in_file, delimiter='\t')
        header = next(reader)
        new_header = [e for e in header]
        index = header.index('Junction.nucleotide.sequence')
        new_header.insert(index+1, 'Generation.probability')
        new_header = '\t'.join(new_header) + '\n'
        outfile = [new_header]
        for line in reader:
            cdr3 = line[header.index('Junction.nucleotide.sequence')]
            p = pgen.compute_nt_CDR3_pgen(cdr3)
            new_line = [e for e in line]
            new_line.insert(index+1, str(p))
            new_line = '\t'.join(new_line) + '\n'
            outfile.append(new_line)
    filename = filename.rstrip('.tsv') + ' (1).tsv'
    with open(filename, 'w') as out_file:
        out_file.writelines(outfile)
    return


def check_thread_status(t, filename):
    coding_pass = False
    while ~coding_pass:
        animate()
        if not t.is_alive():
            break
    print(f'Done processing file \"{filename}\"')
    return


if __name__ == '__main__':
    files = [
        '..\\data_files\\TdTKO\\filtered_data\\filtered_data_exp_TdTKO.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_exp_Normal.tsv',
        # '..\\data_files\\TdTKO\\filtered_data\\filtered_data_gen_TdTKO.tsv',
        # '..\\data_files\\Normal\\filtered_data\\filtered_data_gen_Normal.tsv'
    ]
    processes = []
    pgen = load_olga_mouse_trb_pgen_model()
    for filepath in files:
        file = filepath.split("\\")[-1]
        print(f'Starting to calculate generation probabilities for sequences in file \"{file}\"')
        thread = threading.Thread(target=calculate_pgen, args=(filepath, pgen))
        loading_thread = threading.Thread(target=check_thread_status, args=(thread, file))
        processes.append(loading_thread)
        thread.start()
        loading_thread.start()

    while True in [lt.is_alive() for lt in processes]:
        time.sleep(1)
    print('Done with all files!')
