import csv

v_cdr3 = ['TGCACCTGCAGTGCAGATCTG', 'TGTGCCAGCAGCCAAGATCTT', 'TGTGCCAGCAGCTTAGCGCTA', 'TGTGCCAGCAGCTAAGATCTT',
           'TGTGCCAGCAGCCAAGATCTT', '', '', '', '', '', '', 'TGTGCCAGCTCTCTCGAGA', 'TGTGCCAGCAGTGATGCATC',
           'TGTGCCAGCTCTCTCGAGA', 'TGTGCCAGCGGTGATGCATC', '', 'TGTGCCAGCAGTGATGCATC',
           'TGTGCCAGCAGTTTCTAGAA', 'TGTGCCAGCAGTTTAGCGCTA', 'TGTGCAAGCAGCTTAGATCTA', 'TGTGCTAGCAGTAGAGATCTC',
           '', 'TGTGCCAGCAGTATAGCTAT', 'TGTGGTGCTAGGGATCCC', '', '', 'TGCTCCAGCAGTCAATCGATT', 'TGTGCCAGCAGTCTGTATACA',
           '', 'TGTGCCAGCAGTCTGTCGACA', '', '', 'TGTGCTAGCAGTTTATCGATA', 'TGTAGTTCTAGAGATCTC', 'TGTGCCTGGAGTCTAGAC']
j_cdr3 = ['TTTGCAAACACAGAAGTCTTCTTT', 'TTTGCAAACTCCGACTACACCTTC', 'AGAATTCTGGAAATACGCTCTATTTT',
          'GAAATTTCCAACGAAAGATTATTTTTC', 'GTTATAACAACCAGGCTCCGCTTTTT', '', '', 'GTTATAACTATGCTGAGCAGTTCTTC',
          'TTTGCAAACACCGGGCAGCTCTACTTT', 'CACTAGTGCAGAAACGCTGTATTTT', 'GACTAGTCAAAACACCTTGTACTTT',
          'GGTTAACCAAGACACCCAGTACTTT', '', 'GGAGCTCCTATGAACAGTACTTC']
d_cdr3 = ['TCCCGGGACAGGGGGCGCCC', 'TCCCGGGACTGGGGGGGCGCCC']

with open('..\\..\\data_files\\mus_musculus_vdj_gene_segments.tsv', 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    for line in reader:
        gene_id = line[1]
        gene_seq = line[5]
        cdr3_index = line[3]
        passed = False
        if 'V' in gene_id:
            gene_cdr3 = gene_seq[int(cdr3_index) - 3:]
            print(gene_cdr3)
            for s in v_cdr3:
                cdr3 = s[:-4]
                if gene_cdr3 == cdr3:
                    olga_index = v_cdr3.index(s)
                    print(f'olga_id:{olga_index}\tgene_id:{gene_id}')
                    passed = True
                    break
            if not passed:
                for s in v_cdr3:
                    if gene_cdr3 == s:
                        olga_index = v_cdr3.index(s)
                        print(f'olga_id:{olga_index}\tgene_id:{gene_id}')
                        passed = True
                        break

        elif 'J' in gene_id:
            gene_cdr3 = gene_seq[:int(cdr3_index) + 3]
            print(gene_cdr3)
            for s in j_cdr3:
                cdr3 = s[4:]
                if gene_cdr3 == cdr3:
                    olga_index = j_cdr3.index(s)
                    print(f'olga_id:{olga_index}\t gene_id:{gene_id}')
                    break
        else:
            gene_cdr3 = gene_seq
            print(gene_cdr3)
            for s in d_cdr3:
                cdr3 = s[4:-4]
                if gene_cdr3 == cdr3:
                    olga_index = d_cdr3.index(s)
                    print(f'olga_id:{olga_index}\tgene_id:{gene_id}')
                    break
