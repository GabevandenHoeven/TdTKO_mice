import sys
import time
import difflib
import csv
import re


def animate():
    sys.stdout.write('Loading ')
    for _ in range(4):
        time.sleep(.5)
        sys.stdout.write('.')
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.write('\r')


def check_complementary(nt1: str, nt2: str):
    """Checks if two nucleotides are complementary to each other

    :param nt1: str - A single nucleotide to be checked if it is complementary to the other nucleotide given to
    this function. It can be either DNA or RNA.
    :param nt2: str - The second nucleotide to be checked if it is complementary to the first nucleotide given to
    this function. It can be either DNA or RNA.
    :return:
    """
    if [nt1.upper(), nt2.upper()] in [['A', 'T'], ['T', 'A'], ['G', 'C'], ['C', 'G'], ['A', 'U'], ['U', 'A']]:
        return True
    else:
        return False


def find_p_nucleotides_loop_right(expression, insertion, template, p_nucleotides):
    """

    :param expression:
    :param insertion:
    :param template:
    :param p_nucleotides:
    :return:
    """
    if expression:
        break_ = False
        rev = insertion[::-1]
        for i in range(min(len(rev), len(template))):
            if check_complementary(rev[i], template[i]):
                p_nucleotides += 1
            elif i != 0:
                insertion = insertion[:-i]
                # Removing the found P nucleotides, so they can't be found twice if they match up with the
                # palindromic sequence of the other segment.
                # TODO: if you remove them here there might be a better match to the D segment that gets ignored
                break_ = True
                break
            else:
                break_ = True
        if not break_:
            insertion = ''

    return p_nucleotides, insertion


def find_p_nucleotides_loop_left(expression, insertion, template, p_nucleotides):
    """This function reverses the non-insertion sequence, and then for the length of the shortest sequence between
    the insertion and the template, tries to match every nucleotide as complementary to the one in the other sequence.


    :param expression:
    :param insertion:
    :param template:
    :param p_nucleotides:
    :return:
    """
    if expression:
        break_ = False
        rev = template[::-1]
        for i in range(min(len(insertion), len(rev))):
            if check_complementary(insertion[i], rev[i]):
                p_nucleotides += 1
            else:
                insertion = insertion[i:]
                # Removing the found P nucleotides, so they can't be found twice if they match up with the
                # palindromic sequence of the other segment.
                break_ = True
                break
        if not break_:
            insertion = ''

    return p_nucleotides, insertion


def get_vdj_lengths(input_list: list, RTCR_ref: dict):
    """Code from P.C. de Greef: Accepts a list containing a V- and J-allele and a CDR3 sequence,
    and a dictionary containing reference VDJ alleles with their respective nucleotide sequences.
    This function matches the given alleles to the CDR3, and searches for remnants of the D sequence and likely
    inserted nucleotides. Also checks the non-matching nucleotides for palindromic nucleotides.
    Adapted by Gabe van den Hoeven May 2024.

    :param input_list: list - [v_gene, j_gene, cdr3 nucleotide sequence]
    :param RTCR_ref: dict - dictionary with all the TCR reference sequences, in the function read_rtcr_refs it can be
    specified which organism.
    :returns line_result: list -
    [#nt used of V, matched D, #nt matching D, #nt used of J, #nt deleted of V, #nt deleted of J, #nt left insertions,
    #nt right insertions, #nt left p, #nt right p, last nt of V in the CDR3 sequence,
    first nt of J in the CDR3 sequence]
    """
    V, J, CDR3nt = input_list
    # the specific V and J genes and the whole sequence combined
    if V not in RTCR_ref or J not in RTCR_ref:
        return [25] * 6  # not possible
    germline_V = RTCR_ref[V][0][RTCR_ref[V][1] - 3:]
    # from the allele sequence get the germline sequence, up until the end
    germline_J = RTCR_ref[J][0][:RTCR_ref[J][1] + 3]
    # from the allele sequence get the germline sequence, up until the germline indicator

    match_V = -1
    for i in range(min(len(germline_V), len(CDR3nt))):
        # Match the start of CDR3 to the end of V (germline)
        if germline_V[i] == CDR3nt[i]:
            match_V = i
        else:
            break
    if match_V == -1:
        Vdel = len(germline_V)
        noV_CDR3 = CDR3nt
        Vused = ""
    else:
        Vdel = len(germline_V) - match_V - 1
        noV_CDR3 = CDR3nt[match_V + 1:]
        Vused = germline_V[:match_V + 1]

    match_J = 0
    for i in range(1, min(len(germline_J), len(noV_CDR3)) + 1):
        # Match the end of CDR3 to the start of J (germline)
        if germline_J[-i] == noV_CDR3[-i]:
            match_J = i
        else:
            break
    if match_J == 0:
        Jdel = len(germline_J)
        noVJ_CDR3 = noV_CDR3
        Jused = ""
    else:
        Jdel = len(germline_J) - match_J
        noVJ_CDR3 = noV_CDR3[:-match_J]
        Jused = germline_J[-match_J:]

    D_seq_string = "GGGACAGGGGGC,GGGACTGGGGGGGC".upper()
    # mice seqs from IMGT

    # This is TRB, so check for D
    d = difflib.SequenceMatcher(None, noVJ_CDR3, D_seq_string)
    matchingDlen = max(d.get_matching_blocks(), key=lambda x: x[2])[2]
    d_match = d.find_longest_match()
    Dused = noVJ_CDR3[d_match.a:d_match.a + matchingDlen]
    inslen = len(noVJ_CDR3) - matchingDlen

    # if CDR3nt == 'TGTGCCAGCAGTCCGTCCCGGGACGAAAGATTATTTTTC':
    # if CDR3nt == 'TGTGCAAGCAGCTTAGATGGGACATATGAACAGTACTTC':
    #     print()

    # Checking for palindromic nt, not all insertions are n-nucleotides
    lp_nucleotides, rp_nucleotides, l_ins, r_ins = 0, 0, '', ''
    if inslen != 0:
        p_nt_pairs = []
        ins_pairs = []
        matches = re.finditer(rf'({Dused})', noVJ_CDR3)
        for match in matches:
            lp_nucleotides, rp_nucleotides, l_ins, r_ins = 0, 0, '', ''
            l_ins, r_ins = noVJ_CDR3[:match.start()], noVJ_CDR3[match.end():]
            ins_pairs.append((l_ins, r_ins))

            # Check between D and J
            p_nucleotides1, r_ins1 = find_p_nucleotides_loop_right((r_ins != '' and Jdel == 0), r_ins, Jused, 0)

            # The palindromic nucleotides can also come from the D segment, but only if there are no deletions.
            p_nucleotides1, r_ins1 = find_p_nucleotides_loop_left(
                (D_seq_string.split(',')[0].endswith(Dused) or D_seq_string.split(',')[1].endswith(Dused))
                and r_ins1 != '', r_ins1, Dused, p_nucleotides1)

            # There might be a better match if you check the D first
            p_nucleotides2, r_ins2 = find_p_nucleotides_loop_left(
                (D_seq_string.split(',')[0].endswith(Dused) or D_seq_string.split(',')[1].endswith(Dused))
                and r_ins != '', r_ins, Dused, 0)
            p_nucleotides2, r_ins2 = find_p_nucleotides_loop_right((r_ins2 != '' and Jdel == 0),
                                                                   r_ins2, Jused, p_nucleotides2)
            rp_nucleotides = max(p_nucleotides1, p_nucleotides2)
            if p_nucleotides1 != p_nucleotides2:
                print()
            # Check between V and D
            p_nucleotides1, l_ins1 = find_p_nucleotides_loop_left(l_ins != '' and Vdel == 0, l_ins, Vused, 0)

            # Since the DJ junction is done first, l_ins can contain palindromic nucleotides from the right hand
            # insertion and the J gene. r_ins is retrieved from ins_pairs
            p_nucleotides1, l_ins1 = find_p_nucleotides_loop_right(
                (D_seq_string.split(',')[0].startswith(Dused) or D_seq_string.split(',')[1].startswith(Dused))
                and l_ins1 != '', l_ins1, Dused + ins_pairs[-1][1] + Jused, p_nucleotides1)

            p_nucleotides2, l_ins2 = find_p_nucleotides_loop_right(
                (D_seq_string.split(',')[0].startswith(Dused) or D_seq_string.split(',')[1].startswith(Dused))
                and l_ins != '', l_ins, Dused + ins_pairs[-1][1] + Jused, 0)
            p_nucleotides2, l_ins2 = find_p_nucleotides_loop_left(l_ins2 != '' and Vdel == 0,
                                                                  l_ins2, Vused, p_nucleotides2)
            lp_nucleotides = max(p_nucleotides1, p_nucleotides2)
            if p_nucleotides1 != p_nucleotides2:
                print()
            p_nt_pairs.append((lp_nucleotides, rp_nucleotides))

        sum_p_nt_matches = [sum(i) for i in p_nt_pairs]
        num_max_sum_size = sum_p_nt_matches.count(max(sum_p_nt_matches))
        if num_max_sum_size == 1:
            lp_nucleotides, rp_nucleotides = p_nt_pairs[sum_p_nt_matches.index(max(sum_p_nt_matches))]
            l_ins, r_ins = ins_pairs[sum_p_nt_matches.index(max(sum_p_nt_matches))]
        else:
            # If there are multiple pairs of lp and rp with the same total number of palindromic nucleotides,
            # the pair with the longest single palindromic sequence is chosen.
            # TODO: longest single palindromic sequence may be the best option here anyway. 2-0, 3-0, 4-0 is better
            #  than 1-1, 2-1, 2-2 etc.
            max_p_match = [max(i) for i in p_nt_pairs]
            lp_nucleotides, rp_nucleotides = p_nt_pairs[max_p_match.index(max(max_p_match))]
            l_ins, r_ins = ins_pairs[max_p_match.index(max(max_p_match))]

    line_result = [str(len(Vused)), Dused, str(matchingDlen), str(len(Jused)), str(Vdel), str(Jdel), str(len(l_ins)),
                   str(len(r_ins)), str(lp_nucleotides), str(rp_nucleotides), str(match_V + 1),
                   str(len(Vused) + len(noVJ_CDR3) + 1)]
    return line_result


def read_rtcr_refs():
    """Parses a file with reference alleles from RTCR. Filters for relevant alleles and stores those in a dictionary.
    By P.C. De Greef, 2021
    :returns RTCR_ref: dict - A dictionary containing TCR segment sequences of the beta-chain of Mus Musculus, together
    with the start or stop position of the CDR3 region in the respective segment.
    """
    # read RTCR refs
    # RTCR_ref_filename = "../data_files/immune_receptor_reference.tsv"
    # replace with the gunzipped download of https://github.com/uubram/RTCR/tree/master/rtcr
    RTCR_ref_filename = "..\\data_files\\immune_receptor_reference.tsv"

    RTCR_ref = {}
    with open(RTCR_ref_filename, "r") as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        header = next(reader)
        for row in reader:
            if row[0] == "MusMusculus" and "TR" in row[1]:
                if row[1] not in RTCR_ref:
                    RTCR_ref[row[1]] = [row[-1], int(row[-3])]
                else:
                    print(row[0], "observed twice")
    return RTCR_ref


def reconfigure_header_tdt_normal_data(fn):
    """Used to fix the header of the wild type data files, the header was missing a column making it difficult
    to use it to get values from a specific column.
    :param fn: The path of the file that needed to be edited.

    :return:
    """
    with open(fn, 'r') as file:
        header = next(file)
        header = header.split('\t')
        header.insert(0, 'Column1')
        header = "\t".join(header)
        out = []
        # TODO: readlines?
        for row in file:
            out.append(row)
    print('Reconfigured header.\nWriting to file...')
    with open(fn, 'w') as file:
        file.write(header)
        file.writelines(out)
