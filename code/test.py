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
    complementary = []
    for index in range(len(rev)):
        match rev[index]:
            case 'A':
                complementary.append('T')
            case 'T':
                complementary.append('A')
            case 'C':
                complementary.append('G')
            case 'G':
                complementary.append('C')

    complementary = ''.join(complementary)
    print(complementary)

    n = [1, 2, 3, 4] + ['string']
    print(n)

    sequences = ['a', 'b', 'b', 'c', 'd', 'e', 'd', 'd']
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
    # exp = 'd <= 2'
    # counts = {}
    # for sequence in sequences:
    #     count = sequences.count(sequence)
    #     d = int(d_lengths[sequences.index(sequence)])
    #     if count in range(1, 14) and eval(exp):
    #         mouse = mice[sequences.index(sequence)]
    #         try:
    #             if mouse not in counts[count][2]:
    #                 counts[count][0] += 1
    #                 counts[count][1] += d
    #                 counts[count][2].append(mouse)
    #         except KeyError:
    #             counts.update({count: [1, d, [mouse]]})
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


    def most_frequent(List):
        occurence_count = Counter(List)
        return occurence_count.most_common(1)[0][0]


    List = [2, 1, 2, 2, 1, 3]
    print(most_frequent(List))


def plot():
    plt.figure()
    x = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
    y = [[6.35, 6.30, 6.25, 6.10, 6.05, 5.97, 5.94, 5.88, 5.70, 5.65, 5.50, 5.60, 5.35],
         [6.29, 6.01, 5.94, 5.55, 5.01, 4.95, 4.88, 4.36, 4.23, 3.75]]

    plt.plot(x[0], y[0], label='TdTKO', color='blue')
    plt.plot(x[1], y[1], label='Normal', color='orange')
    plt.scatter(x[1], y[1], color='orange')
    plt.scatter(x[0], y[0], color='blue')
    plt.xticks(numpy.arange(1, 14))
    plt.title('Mean inferred D-length')
    plt.xlabel('Inference')
    plt.ylabel('Inferred D-segment length (nt)')
    plt.legend()

    plt.show()
    plt.close()


plot()

