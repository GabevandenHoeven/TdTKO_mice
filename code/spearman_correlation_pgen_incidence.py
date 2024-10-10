from utils import get_unique_sequences_per_mouse_from_file
from abundance_sequences import get_abundance
import pandas
from scipy.stats import spearmanr
from plots import plot_pgens


if __name__ == '__main__':
    files = [
        '..\\data_files\\WT\\filtered_data\\filtered_data_exp_WT_v2 (1).tsv'
        # '..\\data_files\\WT\\filtered_data\\filtered_data_gen_WT_v2.tsv'
    ]
    mice_per_file = 10
    # mice_per_file = 20

    for file in files:
        max_incidence = mice_per_file
        filtered_data = get_unique_sequences_per_mouse_from_file(file)
        fract_incidence, seq_per_incidence = get_abundance(filtered_data, 'seq[3] >= 0', max_incidence)
        data = []
        data_per_incidence = []
        x = []
        for i in range(len(seq_per_incidence)):
            incidence_group = seq_per_incidence[i]
            for sequence in incidence_group:
                pgen = sequence[10]
                data.append([(i+1) / max_incidence, pgen])
            x.append((i+1) / max_incidence)
            data_per_incidence.append([sequence[10] for sequence in incidence_group])
        newdf = pandas.DataFrame(data, columns=['Incidence', 'Pgen'])
        print(newdf.corr(method='spearman'))
        print(spearmanr(newdf))
        title = 'Generation probabilities per incidence'
        outfile = '..\\img\\Pgen_per_incidence (1).png'
        labels = ['Fraction of incidence', 'Pgen']
        plot_pgens(data_per_incidence, x, title, labels, outfile)
