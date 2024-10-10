import sys
import time
import difflib
import csv
import re
import statistics
import math
from statistics import StatisticsError


def sequence_matcher(a: str, b: str):
    """For every nucleotide in the substring a, this function tries to match the remainder of the substring to
    the super strings: b. The starting indexes in the substring of all unique matches will be returned together with
    the length of the match. If the nucleotide in a is not found in b, no match will be returned.

    :param a: str - The substring.
    :param b: str - The super strings. If there are multiple possibilities they should be separated using ','.
    :return: All unique matches between the remainder of a and any part in b for every nucleotide in a

    Examples:
    a: CGGT
    b: GGGACAGGGGGC,GGGACTGGGGGGGC
    returns: [(0, 1), (1, 1), (1, 2), (2, 1), (3, 1)]

    a: CGGGGT
    b: GGGACAGGGGGC,GGGACTGGGGGGGC
    returns: [(0, 1), (1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (4, 1), (5, 1)]

    a: CGGTAT
    b: GGGACAGGGGGC,GGGACTGGGGGGGC
    returns: [(0, 1), (1, 1), (1, 2), (2, 1), (3, 1), (4, 1), (5, 1)]
    """
    matches = []

    for i in range(len(a)):
        start_indexes = [m.start() for m in re.finditer(a[i], b)]
        match_lengths = [1] * len(start_indexes)

        for j in range(len(start_indexes)):
            ii = i
            k = start_indexes[j]
            coding_pass = True
            while coding_pass:
                try:
                    for m in range(len(match_lengths)):
                        if (i, match_lengths[m]) not in matches:
                            matches.append((i, match_lengths[m]))
                    if a[ii + 1] == b[k + 1]:
                        k += 1
                        ii += 1
                        match_lengths[j] += 1
                    else:
                        coding_pass = False
                except IndexError as e:
                    if e.args[0] == 'string index out of range':
                        coding_pass = False
                        pass
                    else:
                        raise IndexError(e.args)

    return matches


def animate():
    sys.stdout.write('Loading ')
    for _ in range(4):
        time.sleep(.5)
        sys.stdout.write('.')
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.write('\r')


def check_complementary(nt1: str, nt2: str):
    """Checks if two nucleotides are complementary to each other.

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
    """This function tries to match nucleotides from a template to an insertion sequence as being complementary.
    This function reverses the template sequence like for the V segment.

    example:
    V: NNNNNNTGCCGC
    ins: GCGGCANNNNNN

    V becomes CGCCGTNNNNNN

    template    CGCCGTNNNNNN
    insertion   GCGGCANNNNNN

    These sequences match.

    :param expression: bool - This boolean is based on an expression of the starting nucleotides or if there was
    a deletion.
    :param insertion: str - The sequence that is to be compared to the template.
    :param template: str - The nucleotide sequence of the V, potential D, or J segment.
    :param p_nucleotides: int - The number of P nucleotides that were already found.
    :returns p_nucleotides: int - The new number of P nucleotides found.
    :returns insertion: str - The insertion sequence that is left after trying to find P nucleotides.
    """
    if expression:
        break_, i = False, 0
        rev = insertion[::-1]
        for i in range(min(len(rev), len(template))):
            if check_complementary(rev[i], template[i]):
                p_nucleotides += 1
            elif i == 0:
                # Because this function checks for the end of the insertion, if there is no match no nucleotides should
                # be removed. This would happen in the else statement.
                break_ = True
                break
            else:
                insertion = insertion[:-i]
                # Removing the found P nucleotides, so they can't be found twice if they match up with the
                # palindromic sequence of the other segment.
                break_ = True
                break
        if not break_:
            insertion = insertion[:-(i + 1)]

    return p_nucleotides, insertion


def find_p_nucleotides_loop_left(expression, insertion, template, p_nucleotides):
    """This function reverses the non-insertion sequence, and then for the length of the shortest sequence between
    the insertion and the template, tries to match every nucleotide as complementary to the one in the other sequence.

    This function works similar to `find_p_nucleotides_loop_right`


    :param expression: bool - This boolean is based on an expression of the starting nucleotides or if there was
    a deletion.
    :param insertion: str - The sequence that is to be compared to the template.
    :param template: str - The nucleotide sequence of the V, potential D, or J segment.
    :param p_nucleotides: int - The number of P nucleotides that were already found.
    :returns p_nucleotides: int - The new number of P nucleotides found.
    :returns insertion: str - The insertion sequence that is left after trying to find P nucleotides.

    """
    if expression:
        break_, i = False, 0
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
            insertion = insertion[i + 1:]

    return p_nucleotides, insertion


def get_vdj_lengths(input_list: list, RTCR_ref: dict):
    """Code from P.C. de Greef: Accepts a list containing a v- and j-allele and a CDR3 sequence,
    and a dictionary containing reference VDJ alleles with their respective nucleotide sequences.
    This function matches the given alleles to the CDR3, and searches for remnants of the D sequence and likely
    inserted nucleotides. Also checks the non-matching nucleotides for palindromic nucleotides.
    Adapted by Gabe van den Hoeven May 2024.

    :param input_list: list - [v_gene, j_gene, cdr3 nucleotide sequence]
    :param RTCR_ref: dict - dictionary with all the TCR reference sequences, in the function read_rtcr_refs it can be
    specified which organism.
    :returns line_result: list -
    [#nt used of V, matched D, #nt matching D, #nt used of J, #nt deleted of V, #nt deleted of J, #nt left insertions,
    #nt right insertions, #nt left P, #nt right P, last nt of V in the CDR3 sequence,
    first nt of J in the CDR3 sequence]
    """
    v, j, cdr3nt = input_list
    # the specific v and j genes and the whole sequence combined
    if v not in RTCR_ref or j not in RTCR_ref:
        return [25] * 6  # not possible
    germline_v = RTCR_ref[v][0][RTCR_ref[v][1] - 3:]
    # from the allele sequence get the germline sequence, up until the end
    germline_j = RTCR_ref[j][0][:RTCR_ref[j][1] + 3]
    # from the allele sequence get the germline sequence, up until the germline indicator

    match_v = -1
    for i in range(min(len(germline_v), len(cdr3nt))):
        # Match the start of CDR3 to the end of v (germline)
        if germline_v[i] == cdr3nt[i]:
            match_v = i
        else:
            break
    if match_v == -1:
        v_del = len(germline_v)
        no_v_cdr3 = cdr3nt
        v_used = ""
    else:
        v_del = len(germline_v) - match_v - 1
        no_v_cdr3 = cdr3nt[match_v + 1:]
        v_used = germline_v[:match_v + 1]

    match_j = 0
    for i in range(1, min(len(germline_j), len(no_v_cdr3)) + 1):
        # Match the end of CDR3 to the start of j (germline)
        if germline_j[-i] == no_v_cdr3[-i]:
            match_j = i
        else:
            break
    if match_j == 0:
        j_del = len(germline_j)
        no_vj_cdr3 = no_v_cdr3
        j_used = ""
    else:
        j_del = len(germline_j) - match_j
        no_vj_cdr3 = no_v_cdr3[:-match_j]
        j_used = germline_j[-match_j:]

    d_seq_string = "GGGACAGGGGGC,GGGACTGGGGGGGC".upper()
    # mice seqs from IMGT

    # This is TRB, so check for D
    matching_d_len, lp_nucleotides, rp_nucleotides, d_used, l_ins, r_ins = 0, 0, 0, '', '', ''

    possible_ds = sequence_matcher(no_vj_cdr3, d_seq_string)
    d_p_matches = []

    for d_match in possible_ds:
        matching_d_len = d_match[1]
        d_used = no_vj_cdr3[d_match[0]:d_match[0] + matching_d_len]
        ins_len = len(no_vj_cdr3) - matching_d_len
        lp_nucleotides, rp_nucleotides, l_ins, r_ins = 0, 0, '', ''
        # Checking for palindromic nt, not all insertions are n-nucleotides
        if ins_len != 0:
            l_ins, r_ins = no_vj_cdr3[:d_match[0]], no_vj_cdr3[d_match[0] + d_match[1]:]

            # Check between D and j
            p_nucleotides1, r_ins1 = find_p_nucleotides_loop_right((r_ins != '' and j_del == 0), r_ins, j_used, 0)

            # The palindromic nucleotides can also come from the D segment, but only if there are no deletions.
            p_nucleotides1, r_ins1 = find_p_nucleotides_loop_left(
                (d_seq_string.split(',')[0].endswith(d_used) or d_seq_string.split(',')[1].endswith(d_used))
                and r_ins1 != '', r_ins1, d_used, p_nucleotides1)

            # There might be a better match if you check the D first
            p_nucleotides2, r_ins2 = find_p_nucleotides_loop_left(
                (d_seq_string.split(',')[0].endswith(d_used) or d_seq_string.split(',')[1].endswith(d_used))
                and r_ins != '', r_ins, d_used, 0)
            p_nucleotides2, r_ins2 = find_p_nucleotides_loop_right((r_ins2 != '' and j_del == 0),
                                                                   r_ins2, j_used, p_nucleotides2)
            rp_nucleotides = max(p_nucleotides1, p_nucleotides2)

            # Check between v and D
            p_nucleotides1, l_ins1 = find_p_nucleotides_loop_left(l_ins != '' and v_del == 0, l_ins, v_used, 0)

            # Since the DJ junction is done first, l_ins can contain palindromic nucleotides from the right-hand
            # insertion and the j gene. r_ins is retrieved from ins_pairs
            p_nucleotides1, l_ins1 = find_p_nucleotides_loop_right(
                (d_seq_string.split(',')[0].startswith(d_used) or d_seq_string.split(',')[1].startswith(d_used))
                and l_ins1 != '', l_ins1, d_used + r_ins + j_used, p_nucleotides1)

            p_nucleotides2, l_ins2 = find_p_nucleotides_loop_right(
                (d_seq_string.split(',')[0].startswith(d_used) or d_seq_string.split(',')[1].startswith(d_used))
                and l_ins != '', l_ins, d_used + r_ins + j_used, 0)
            p_nucleotides2, l_ins2 = find_p_nucleotides_loop_left(l_ins2 != '' and v_del == 0,
                                                                  l_ins2, v_used, p_nucleotides2)
            lp_nucleotides = max(p_nucleotides1, p_nucleotides2)

        d_p_matches.append((matching_d_len, lp_nucleotides, rp_nucleotides, d_used, l_ins, r_ins))

    try:
        # If the no_vj_cdr3 completely matches a D segment that match is chosen.
        max_d_p_nt = max([(match[0] + match[1] + match[2]) for match in d_p_matches])
        max_d_p_nt_list = [match for match in d_p_matches if (match[0] + match[1] + match[2]) == max_d_p_nt]
        best_match = max(max_d_p_nt_list, key=lambda x: x[0])
        matching_d_len, lp_nucleotides, rp_nucleotides, d_used, l_ins, r_ins = best_match
    except ValueError as e:
        if e.args[0] == 'max() arg is an empty sequence':
            pass
    line_result = [str(len(v_used)), d_used, str(matching_d_len), str(len(j_used)), str(v_del), str(j_del),
                   str(len(l_ins)), str(len(r_ins)), str(lp_nucleotides), str(rp_nucleotides), str(match_v + 1),
                   str(len(v_used) + len(no_vj_cdr3) + 1)]
    return line_result


def get_vdj_lengths_old_method(input_list: list, RTCR_ref: dict):
    """Code from P.C. de Greef: Accepts a list containing a V- and J-allele and a CDR3 sequence,
    and a dictionary containing reference VDJ alleles with their respective nucleotide sequences.
    This function matches the given alleles to the CDR3, and searches for remnants of the D sequence and likely
    inserted nucleotides.
    Adapted by Gabe van den Hoeven May 2024.

    :param input_list: list - [v_gene, j_gene, cdr3 nucleotide sequence]
    :param RTCR_ref: dict - dictionary with all the TCR reference sequences, in the function read_rtcr_refs it can be
    specified which organism.
    :returns line_result: list -
    [#nt used of V, matched D, #nt matching D, #nt used of J, #nt deleted of V, #nt deleted of J, last nt of V in
    the CDR3 sequence, first nt of J in the CDR3 sequence]
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
    matchingDlen, Dused = 0, ''

    d_gene = difflib.SequenceMatcher(None, noVJ_CDR3, D_seq_string)
    matchingDlen = max(d_gene.get_matching_blocks(), key=lambda x: x[2])[2]
    d_match = max(d_gene.get_matching_blocks(), key=lambda x: x[2])
    Dused = noVJ_CDR3[d_match.a: d_match.a + matchingDlen]
    l_ins, r_ins = noVJ_CDR3[:d_match.a], noVJ_CDR3[d_match.a + matchingDlen:]

    line_result = [str(len(Vused)), Dused, str(matchingDlen), str(len(Jused)), str(Vdel), str(Jdel),
                   str(len(l_ins) + len(r_ins)), str(match_V + 1), str(len(Vused) + len(noVJ_CDR3) + 1)]
    return line_result


def read_rtcr_refs():
    """Parses a data with reference alleles from RTCR. Filters for relevant alleles and stores those in a dictionary.
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


def reconfigure_header_tdt_wt_data(fn):
    """Used to fix the header of the wild type data files, the header was missing a column making it difficult
    to use it to get values from a specific column.
    :param fn: The path of the data that needed to be edited.

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
    print('Reconfigured header.\nWriting to data...')
    with open(fn, 'w') as file:
        file.write(header)
        file.writelines(out)


def calculate_confidence_intervals(values, z=1.96):
    """This function calculates the confidence intervals of a given set of values.
    Returns the mean and confidence intervals.

    :param values: Iterable obj - The raw values for which to calculate confidence intervals.
    :param z: float - The Z-score for calculating the confidence intervals. Default = 1.96
    :return:
    """
    mean = statistics.mean(values)
    try:
        stdev = statistics.stdev(values)
    except StatisticsError:
        stdev = 0
    confidence_interval = z * (stdev / math.sqrt(len(values)))
    return mean, confidence_interval


def get_unique_sequences_from_file(filename, delim='\t'):
    """Checks a data data for the CDR3 sequences in it to be unique. The first line of a unique sequence in the data,
    gets added to the output list.
    For example:
    Column1	Mouse	Number.of.reads	V.length.used	D.used	D.length.used	J.length.used	V.length.deleted
    J.length.deleted	Left.insertion.length	Right.insertion.length	Left.palindromic	Right.palindromic
    V.J.distance	V.gene	J.gene	V.gene.end.position	J.gene.start.position	Min.Phred
    Junction.nucleotide.sequence	phenotype	strain
    3	M28	6682	9		0	21	6	0	0	0	0	0	0	TRBV12-2*01	TRBJ2-3*01	9	10	40	TGTGCCAGCAGTGCAGAAACGCTGTATTTT	CD4+	TdT-/-
    6	M28	5623	9		0	21	6	0	0	0	0	0	0	TRBV12-2*01	TRBJ2-3*01	9	10	40	TGTGCCAGCAGTGCAGAAACGCTGTATTTT	CD4+	TdT-/-
    8	M28	2761	17		0	16	0	4	0	0	0	0	0	TRBV16*01	TRBJ1-1*01	17	18	40	TGTGCAAGCAGCTTAGACACAGAAGTCTTCTTT	CD4+	TdT-/-
    24	M29	1966	9		0	21	6	0	0	0	0	0	0	TRBV12-2*01	TRBJ2-3*01	9	10	40	TGTGCCAGCAGTGCAGAAACGCTGTATTTT	CD4+	TdT-/-

    Here, line 3, 8 and 24 get added because they are unique sequences or belong to a different mouse,
    but line 6 is the same as line 3 in sequence and mouse, so it is disregarded.

    :param filename: str - The name of the data, should be a path if it is not within the current working directory
    :param delim: str - The delimiter of the datafile, which character is used to separate values.
    :return:
    """
    new_lines, sequences_check = [], {}
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=delim)
        header = next(reader)
        new_lines.append(header)
        for line in reader:
            try:
                sequences_check[
                    line[header.index('V.gene')] +
                    line[header.index('Junction.nucleotide.sequence')] +
                    line[header.index('J.gene')]]
            except KeyError:
                new_lines.append(line)
                sequences_check.update({line[header.index('V.gene')] +
                                        line[header.index('Junction.nucleotide.sequence')] +
                                        line[header.index('J.gene')]: True})
    return new_lines


def get_unique_sequences_per_mouse_from_file(filename, delim='\t'):
    """This function is a variation on the ``get_unique_sequences_from_file`` function.
    This function allows for the same sequence to be added if the sequence is found in a different mouse.

    :param filename: str - The name of the data, should be a path if it is not within the current working directory
    :param delim: str - The delimiter of the datafile, which character is used to separate values.
    :return:
    """
    new_lines, sequences_check, dups = [], {}, 0
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter=delim)
        header = next(reader)
        new_lines.append(header)
        for line in reader:
            try:
                if sequences_check[
                                    line[header.index('Mouse')] +
                                    line[header.index('V.gene')] +
                                    line[header.index('Junction.nucleotide.sequence')] +
                                    line[header.index('J.gene')]]:
                    dups += 1
            except KeyError:
                new_lines.append(line)
                sequences_check.update({line[header.index('Mouse')] +
                                        line[header.index('V.gene')] +
                                        line[header.index('Junction.nucleotide.sequence')] +
                                        line[header.index('J.gene')]: True})
    return new_lines


def flatten_list_unique(obj):
    """Flattens a 2-dimensional iterable object into a list.

    :param obj: a 2-dimensional iterable object like a nested list.
    :return:
    """
    def flatten(x):
        if isinstance(x, (list, tuple, set, range)):
            for x_ in x:
                yield from flatten_list_unique(x_)
        else:
            yield x
    return list(flatten(obj))
