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
