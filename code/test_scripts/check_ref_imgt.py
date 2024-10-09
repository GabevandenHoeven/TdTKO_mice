import csv
import re
from parse_imgt import parse_imgt_webpage


def read_rtcr_refs():
    """Parses a data with reference alleles from RTCR. Filters for relevant alleles and stores those in a dictionary.
    By P.C. De Greef, 2021
    :returns RTCR_ref: dict - A dictionary containing TCR segment sequences of the beta-chain of Mus Musculus, together
    with the start or stop position of the CDR3 region in the respective segment.
    """
    # read RTCR refs
    # RTCR_ref_filename = "../data_files/immune_receptor_reference.tsv"
    # replace with the gunzipped download of https://github.com/uubram/RTCR/tree/master/rtcr
    RTCR_ref_filename = "..\\..\\data_files\\immune_receptor_reference.tsv"

    RTCR_ref = {}
    with open(RTCR_ref_filename, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        header = next(reader)
        for row in reader:
            if row[0] == "MusMusculus" and re.search("TRB[VDJ]{1}", row[1]):
                if row[1] not in RTCR_ref:
                    RTCR_ref[row[1]] = [row[-1], int(row[-3])]
                else:
                    print(row[0], "observed twice")
    return RTCR_ref


def get_vdj_segments_from_olga_param_file(filename):
    vdj_segments = {}
    with open(filename, 'r') as infile:
        for line in infile:
            if line.startswith("#GeneChoice"):
                read = True
            elif line.startswith('%') and read:
                segment, sequence = line.lstrip('%').split(';')[0], line.split(';')[1]
                vdj_segments.update({segment: sequence})
            else:
                read = False
    return vdj_segments


if __name__ == '__main__':
    tdt_alleles = ('TRBV1*01', 'TRBV2*01', 'TRBV3*01', 'TRBV4*01', 'TRBV5*01', 'TRBV12-1*01', 'TRBV12-2*01',
                   'TRBV13-1*01', 'TRBV13-2*01', 'TRBV13-3*01', 'TRBV14*01', 'TRBV15*01', 'TRBV16*01', 'TRBV17*01',
                   'TRBV19*01', 'TRBV20*01', 'TRBV23*01', 'TRBV24*01', 'TRBV26*01', 'TRBV29*01', 'TRBV30*01',
                   'TRBV31*01',

                   'TRBJ1-1*01', 'TRBJ1-2*01', 'TRBJ1-3*01', 'TRBJ1-4*01', 'TRBJ1-5*01', 'TRBJ2-1*01', 'TRBJ2-2*01',
                   'TRBJ2-3*01', 'TRBJ2-4*01', 'TRBJ2-5*01', 'TRBJ2-7*01')
    print(f"length list of alleles in data: {len(tdt_alleles)}")

    all_vdj_alleles_olga = get_vdj_segments_from_olga_param_file(
        '..\\models\\Mus+musculus\\TRB\\models\\model_params_no_err.txt')
    print(f"length list of alleles on OLGA: {len(all_vdj_alleles_olga)}")
    # TODO: THESE MIGHT NOT BE THE EXACT CORRECT NAMES EG: TRBV24*01 / TRBV24*02

    refs = read_rtcr_refs()
    print(f"length list of all TRB[VDJ] alleles in ref data: {len(refs.keys())}")

    results = parse_imgt_webpage(url='https://www.imgt.org/ligmdb/view.action?id=IMGT000132')
    print(f"length list of all TRB[VDJ] alleles on IMGT: {len(results.keys())}")

    print(f"Alleles in data missing in peter reference data: {[e for e in tdt_alleles if e not in refs.keys()]}")

    print(f'Alleles from OLGA missing in peter reference data: '
          f'{[e for e in all_vdj_alleles_olga if e not in refs.keys()]}')

    imgt_alleles = sorted(list(results.keys()))
    olga_alleles = sorted(list(all_vdj_alleles_olga.keys()))
    print(f"Allele symbols on IMGT and OLGA match: {imgt_alleles == olga_alleles}")
    olga_diff_segments = [x for x in olga_alleles if x not in set(imgt_alleles)]
    print(olga_diff_segments, 'Allele symbols different in OLGA')
    imgt_diff_segments = [x for x in imgt_alleles if x not in set(olga_alleles)]
    print(imgt_diff_segments, 'Allele symbols different in IMGT')

    print()
    print(f"Matching OLGA and IMGT sequence: TRBV13-1*01/TRBV13-1*02\n"
          f"{all_vdj_alleles_olga['TRBV13-1*01'].upper() == results['TRBV13-1*02'].upper()}")
    print(f"Matching OLGA and IMGT sequence: TRBV24*01/TRBV24*02\n"
          f"{all_vdj_alleles_olga['TRBV24*01'].upper() == results['TRBV24*02'].upper()}")
    print(f"Matching OLGA and IMGT sequence: TRBJ1-4*01/TRBJ1-4*02\n"
          f"{all_vdj_alleles_olga['TRBJ1-4*01'].upper() == results['TRBJ1-4*02'].upper()}")
    print()

    differences = []
    for e in imgt_alleles:
        if e not in olga_diff_segments and e not in imgt_diff_segments:
            print(f'Matching OLGA and IMGT sequence: {e}\n{all_vdj_alleles_olga[e].upper() == results[e].upper()}')
            if all_vdj_alleles_olga[e].upper() != results[e].upper():
                differences.append(e)
    print('\nThe alleles that do not match:', differences)

    differences = []
    print()
    for e in refs.keys():
        print(f"Matching OLGA and reference data sequence: {e}\n"
              f"{all_vdj_alleles_olga[e].upper() == refs[e][0].upper()}")
        if all_vdj_alleles_olga[e].upper() != refs[e][0].upper():
            differences.append(e)
    print('\nThe alleles that do not match:', differences)

    differences = []
    print()
    for e in refs.keys():
        if e not in olga_diff_segments:
            print(f"Matching IMGT and reference data sequence: {e}\n"
                  f"{results[e].upper() == refs[e][0].upper()}")
            if results[e].upper() != refs[e][0].upper():
                differences.append(e)
    print('\nThe alleles that do not match:', differences)

    print()
    print(f"Matching IMGT and reference sequence: TRBV13-1*01/TRBV13-1*02\n"
          f"{refs['TRBV13-1*01'][0].upper() == results['TRBV13-1*02'].upper()}\n"
          f"{refs['TRBV13-1*01'][0].upper()}\n"
          f"{results['TRBV13-1*02'].upper()}")
    print(f"Matching IMGT and reference sequence: TRBV24*01/TRBV24*02\n"
          f"{refs['TRBV24*01'][0].upper() == results['TRBV24*02'].upper()}\n"
          f"{refs['TRBV24*01'][0].upper()}\n"
          f"{results['TRBV24*02'].upper()}")
    print(f"Matching IMGT and reference sequence: TRBJ1-4*01/TRBJ1-4*02\n"
          f"{refs['TRBJ1-4*01'][0].upper() == results['TRBJ1-4*02'].upper()}\n"
          f"{refs['TRBJ1-4*01'][0].upper()}\n"
          f"{results['TRBJ1-4*02'].upper()}")
