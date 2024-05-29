from collections import Counter
import matplotlib.pyplot as plt
import numpy

if __name__ == '__main__':
    seq = 'TTACGTGATGGGTGGAGTCACATTTCTCAGATCCTCTAGGACGGTGAGTCGTGTCCCTGGTCCGAAGAACTGCTCAGCATACCAGTCCGATCCAGCCTCATTAGAATGTAGTCCCAGACATGAGAGAGCCTCTTCCTCAAGGAGGAACCAT'
    rev = seq[::-1]
    print(rev)
    # seq1 = 'AACTACGTAGGGTGGAGTCACATTTCTCAGATCCTCTACAACAATGAGCCGGCTTCCTTCTCCAAAATAGAGCGTATTTCCACTGCTGGCACAGAGAAAAACGGCCATCTCGTTCTTCTACA'
    # seq2 = 'TTGATGCATCCCACCTCAGTGTAAAGAGTCTAGGAGATGTTGTTACTCGGCCGAAGGAAGAGGTTTTATCTCGCATAAAGGTGACGACCGTGTCTCTTTTTGCCGGTAGAGCAAGAAGATGT'
    #
    # for index in range(len(seq1)):
    #     if seq1[index] == 'A' and seq2[index] == 'T':
    #         print('True')
    #     elif seq1[index] == 'T' and seq2[index] == 'A':
    #         print("True")
    #     elif seq1[index] == 'C' and seq2[index] == 'G':
    #         print('True')
    #     elif seq1[index] == 'G' and seq2[index] == 'C':
    #         print('True')
    #     else:
    #         print('False')
    # complementary = []
    # for index in range(len(rev)):
    #     match rev[index]:
    #         case 'A':
    #             complementary.append('T')
    #         case 'T':
    #             complementary.append('A')
    #         case 'C':
    #             complementary.append('G')
    #         case 'G':
    #             complementary.append('C')
    #
    # complementary = ''.join(complementary)
    # print(complementary)
    #
    # n = [1, 2, 3, 4] + ['string']
    # print(n)

    sequences = ['a', 'b', 'b', 'c', 'd_gene', 'e', 'd_gene', 'd_gene']
    d_lengths = ['1', '4', '4', '3', '1', '7', '2', '1']
    mice = ['M1', 'M1', 'M2', 'M3', 'M4', 'M1', 'M2', 'M3']
    incidences = []
    mean_d_lengths = []
    # for i in range(1, 14):
    #     count = 0
    #     mean_d = 0
    #     for sequence in sequences:
    #         if sequences.count(sequence) == i:
    #             count += 1
    #             mean_d += int(d_lengths[sequences.index(sequence)])
    #     incidences.append(count)
    #     try:
    #         mean_d_lengths.append(mean_d / count)
    #     except ZeroDivisionError:
    #         mean_d_lengths.append(0)
    #     # n_seq_per_incidence = sum(1 for c in Counter(sequences).values() if c == i)
    #     # incidences.append(n_seq_per_incidence)
    # print(incidences)
    # print(mean_d_lengths)
    # exp = 'd_gene <= 2'
    # counts = {}
    # for sequence in sequences:
    #     count = sequences.count(sequence)
    #     d_gene = int(d_lengths[sequences.index(sequence)])
    #     if count in range(1, 14) and eval(exp):
    #         mouse = mice[sequences.index(sequence)]
    #         try:
    #             if mouse not in counts[count][2]:
    #                 counts[count][0] += 1
    #                 counts[count][1] += d_gene
    #                 counts[count][2].append(mouse)
    #         except KeyError:
    #             counts.update({count: [1, d_gene, [mouse]]})
    # for c in sorted(counts.keys()):
    #     incidences.append(c)
    #     mean_d_lengths.append(counts[c][1] / c)
    #     try:
    #         mean_d_lengths.append(mean_d / count)
    #     except ZeroDivisionError:
    #         mean_d_lengths.append(0)
    #     # n_seq_per_incidence = sum(1 for c in Counter(sequences).values() if c == i)
    #     # incidences.append(n_seq_per_incidence)
    # print(incidences)
    # print(mean_d_lengths)
    from collections import Counter


    # def most_frequent(List):
    #     occurence_count = Counter(List)
    #     return occurence_count.most_common(1)[0][0]
    #
    #
    # List = [2, 1, 2, 2, 1, 3]
    # print(most_frequent(List))


# def plot():
#     plt.figure()
#     x = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
#          [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
#     y = [[6.35, 6.30, 6.25, 6.10, 6.05, 5.97, 5.94, 5.88, 5.70, 5.65, 5.50, 5.60, 5.35],
#          [6.29, 6.01, 5.94, 5.55, 5.01, 4.95, 4.88, 4.36, 4.23, 3.75]]
#
#     plt.plot(x[0], y[0], label='TdTKO', color='blue')
#     plt.plot(x[1], y[1], label='Normal', color='orange')
#     plt.scatter(x[1], y[1], color='orange')
#     plt.scatter(x[0], y[0], color='blue')
#     plt.xticks(numpy.arange(1, 14))
#     plt.title('Mean inferred D-length')
#     plt.xlabel('Inference')
#     plt.ylabel('Inferred D-segment length (nt)')
#     plt.legend()
#
#     plt.show()
#     plt.close()
#
#
# plot()

seqs = {
    'agagada': [1, [0], ["m28"], [0], [0]],
    'afagaga': [1, [0], ["m28"], [0], [0]],
    'afasgfgaa': [1, [0], ["m28"], [0], [0]],
    'agdadga': [1, [0], ["m28"], [0], [0]],
    'addasgdav': [1, [0], ["m28"], [0], [0]],
    'adgfsdghsys': [1, [0], ["m28"], [0], [0]],
    'dgsfgasysa': [1, [0], ["m28"], [0], [0]],
    'gshagda': [2, [0, 0], ["m28", "m30"], [0, 0], [0, 0]],
    'gshasga': [2, [0, 0], ["m28", "m30"], [0, 0], [0, 0]],
    'gshassdgagga': [2, [0, 0], ["m28", "m30"], [0, 0], [0, 0]],
    'gsadgaadga': [2, [0, 0], ["m28", "m30"], [0, 0], [0, 0]],
    'aghrgada': [1, [0], ["m28"], [0], [0]],
    'ahfshgagba': [1, [0], ["m28"], [0], [0]],
    'agsdgavv': [4, [0, 0, 0, 0], ["m28", "m30", "m29", "m31"], [0, 0, 0, 0], [0, 0, 0, 0]],
    'agsdgv': [3, [0, 0, 0], ["m28", "m30", "m29"], [0, 0, 0], [0, 0, 0]],
    'agsv': [1, [0], ["m28"], [0], [0]],
    'agsdvv': [1, [0], ["m28"], [0], [0]],
    'aavv': [1, [0], ["m28"], [0], [0]],
    'agsdfgsvv': [1, [0], ["m28"], [0], [0]],
}
# for i in range(1, 14):
#     n_seq_per_incidence = sum(1 for c in seqs.values() if c[0] == i)
#     print(f'{i}: {n_seq_per_incidence}')

p_nt_pairs = [(0, 1), (0, 0), (1, 1), (0, 2), (0, 2)]
ins_pairs = [(0, 1), (0, 0), (1, 1), (0, 2), (1, 2)]
sum_p_nt = [sum(i) for i in p_nt_pairs]
m = max([sum(i) for i in p_nt_pairs])
max_p_nt_list = [p for p in p_nt_pairs if sum(p) == m]
max_p = [max(p) for p in max_p_nt_list]
lp_nucleotides, rp_nucleotides = max_p_nt_list[max_p.index(max(max_p))]
l_ins, r_ins = ins_pairs[p_nt_pairs.index((lp_nucleotides, rp_nucleotides))]
print()
