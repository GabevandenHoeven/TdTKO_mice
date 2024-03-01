from code_pc_degreef import read_rtcr_refs
from parse_imgt import parse_imgt_webpage

tdt_alleles = {'TRBJ1-2*01', 'TRBV12-1*01', 'TRBJ1-5*01', 'TRBV14*01', 'TRBV13-1*01', 'TRBJ1-4*01', 'TRBV30*01', 'TRBV4*01', 'TRBV19*01', 'TRBJ2-5*01', 'TRBV23*01', 'TRBV24*01', 'TRBJ2-3*01', 'TRBJ2-2*01', 'TRBV5*01', 'TRBV20*01', 'TRBV1*01', 'TRBJ2-7*01', 'TRBV13-2*01', 'TRBV13-3*01', 'TRBJ1-3*01', 'TRBV12-2*01', 'TRBJ2-1*01', 'TRBJ2-4*01', 'TRBV17*01', 'TRBV2*01', 'TRBV16*01', 'TRBV3*01', 'TRBJ1-1*01', 'TRBV31*01', 'TRBV29*01', 'TRBV15*01', 'TRBV26*01'}
print(f"length list of alleles in data: {len(tdt_alleles)}\n")

print("Peter reference file missing:\n")
refs = read_rtcr_refs()
for e in tdt_alleles:
    if e not in refs.keys():
        print(e)


print("IMGT webpage missing:\n")
results = parse_imgt_webpage(url='https://www.imgt.org/ligmdb/view.action?id=IMGT000132')
for e in tdt_alleles:
    if e not in results.keys():
        print(e)
        print(refs[e][0])

print("Check if reference file sequences match IMGT:\n")
for e in tdt_alleles:
    ref_seq = refs[e][0]
    try:
        imgt_seq = results[e].upper()
        if ref_seq != imgt_seq:
            print(e)
    except KeyError:
        pass

ref_seq = "TTTCCAACGAAAGATTATTTTTCGGTCATGGAACCAAGCTGTCTGTCTTGG"
imgt_seq = "tttccaacgaaagattatttttcggtcatggaaccaagctgtctgtcttgg".upper()
if ref_seq != imgt_seq:
    print("TRBJ1-4*01")

ref_seq = "GAGGCTGCAGTCACCCAAAGCCCTAGAAACAAGGTGACAGTAACAGGAGGAAACGTGACATTGAGCTGTCGCCAGACTAATAGCCACAACTACATGTACTGGTATCGGCAGGACACTGGGCATGGGCTGAGGCTGATCCATTACTCATATGGTGCTGGCAACCTTCGAATAGGAGATGTCCCTGATGGGTACAAGGCCACCAGAACAACGCAAGAAGACTTCTTCCTCCTGCTGGAATTGGCTTCTCCCTCTCAGACATCTTTGTACTTCTGTGCCAGCAGTGATG"
imgt_seq = "gaggctgcagtcacccaaagccctagaaacaaggtgacagtaacaggaggaaacgtgacattgagctgtcgccagactaatagccacaactacatgtactggtatcggcaggacactgggcatgggctgaggctgatccattactcatatggtgctggcaaccttcgaataggagatgtccctgatgggtacaaggccaccagaacaacgcaagaagacttcttcctcctgctggaattggcttctccctctcagacatctttgtacttctgtgccagcagtgatg".upper()
if ref_seq != imgt_seq:
    print("TRBV13-1*01")

ref_seq = "GTTGCTGGAGTAACCCAGACTCCACGATACCTGGTCAAAGAGAAAGGACAGAAAGCACACATGAGCTGTAGTCCTGAAAAAGGGCACACTGCCTTTTACTGGTATCAACAGAACCAGAAACAAGAACTTACATTTTTGATTAGCTTTCGAAATGAAGAAATTATGGAACAAACAGACTTGGTCAAGAAGAGATTCTCAGCTAAGTGTTCCTCGAACTCACGCTGCATCCTGGAAATCCTATCCTCTGAAGAAGACGACTCAGCACTGTACCTCTGTGCCAGCAGTCTGTA"
imgt_seq = "gttgctggagtaacccagactccacgatacctggtcaaagagaaaggacagaaagcacacatgagctgtagtcctgaaaaagggcacactgccttttactggtatcaacagaaccagaaacaagaacttacatttttgattagctttcgaaatgaagaaattatggaacaaacagacttggtcaagaagagattctcagctaagtgttcctcgaactcacgctgcatcctggaaatcctatcctctgaagaagacgactcagcactgtacctctgtgccagcagtctgta".upper()
if ref_seq != imgt_seq:
    print("TRBV24*01")