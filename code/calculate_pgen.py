import threading
import time
import subprocess
from utils import animate
import csv


def calculate_pgen(filename, delim):
    """This function calculates generation probabilities for a datafile containing CDR3 nucleotide sequences and their
    called V and J gene segments using OLGA. OLGA is required to be installed first.
    This function then inserts these generation probabilities after the nucleotide sequences, and writes
    this to a new file together with the other information in the datafile.

    :param filename: str - A file path to a datafile. This datafile should have a header with, for the nucleotide
    sequence, V and J gene segment columns, the following header names: 'Junction.nucleotide.sequence', 'V.gene',
    'J.gene'.
    :param delim: str - The delimiter of the datafile.
    :return:
    """
    with open(filename, 'r') as in_file:
        reader = csv.reader(in_file, delimiter=delim)
        header = next(reader)
        seq_index, v_index, j_index = str(header.index('Junction.nucleotide.sequence')), \
            str(header.index('V.gene')), str(header.index('J.gene'))

    command = ['olga-compute_pgen', '-i', filename, '--lines_to_skip', '1', '--mouseTRB', '-o', 'pgen_output.tsv',
               '--seq_in', seq_index, '--v_in', v_index, '--j_in', j_index]
    subprocess.run(command)

    with open('pgen_output.tsv', 'r') as pgen_file:
        pgen_reader = csv.reader(pgen_file, delimiter='\t')
        generation_probabilities = []
        for line in pgen_reader:
            generation_probabilities.append(line[1])
            
    with open(filename, 'r') as in_file:
        reader = csv.reader(in_file, delimiter=delim)
        header = next(reader)
        new_header = [e for e in header]
        index = header.index('Junction.nucleotide.sequence')
        new_header.insert(index + 1, 'Generation.probability')
        new_header = '\t'.join(new_header) + '\n'
        outfile = [new_header]
        for line in reader:
            try:
                pgen = generation_probabilities[reader.line_num - 2]
            except IndexError:
                break
            line.insert(index + 1, pgen)
            outfile.append('\t'.join(line) + '\n')
    filename = filename.rstrip('.tsv') + '_with_pgens.tsv'
    with open(filename, 'w') as out_file:
        out_file.writelines(outfile)
    return


def check_thread_status(t, filename):
    """A small function to run an animation in parallel while the generation probabilities are being calculated.
    This is to make sure the script has not stopped, since the pgen calculations can take a long time depending on the
    input datafile.
    :param t: Thread - A Thread object that is running the function that calculates generation probabilities
    :param filename: str - File name of the datafile for which Pgens are being calculated.
    """
    coding_pass = False
    while ~coding_pass:
        animate()
        if not t.is_alive():
            break
    print(f'Done processing file \"{filename}\"')
    return


if __name__ == '__main__':
    files = [
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2.tsv',
    ]
    processes = []
    for filepath in files:
        file = filepath.split("\\")[-1]
        print(f'Starting to calculate generation probabilities for sequences in file \"{file}\"')
        thread = threading.Thread(target=calculate_pgen, args=(filepath, '\t'))
        loading_thread = threading.Thread(target=check_thread_status, args=(thread, file))
        processes.append(loading_thread)
        thread.start()
        loading_thread.start()

    while True in [lt.is_alive() for lt in processes]:
        time.sleep(1)
    print('Done with all files!')
