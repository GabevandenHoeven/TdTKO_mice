import urllib.request
import re


def parse_imgt_webpage(url):
    """
    """
    # urllib.request.urlretrieve(url, 'C:\\Users\\gabev\\PycharmProjects\\MRP\\data_files\\temp.html')
    with open('..\\..\\data_files\\temp.html', 'r') as file:
        ls = ['V-REGION', 'D-REGION', 'J-REGION']
        results = {}
        b = False
        for line in file:
            if any(subs in line for subs in ls):
                b = True
            elif b:
                if re.search('IMGT_allele', line):
                    allele = re.search('(TRB[VDJ][0-9\-\*]+)', next(file)).group()
                elif re.search('>[atgc]+<', line):
                    sequence = re.search('>[atgc]+<', line).group().lstrip('>').rstrip('<')
                    results.update({allele: sequence})
                    b = False

    with open('parsed.tsv', 'w') as file:
        file.write("Allele\tN_sequence\n")
        for allele in sorted(results.keys()):
            file.write(f"{allele}\t{results[allele]}\n")
    return results


if __name__ == '__main__':
    print(re.search('TRB[VDJ][0-9\-\*]+', 'TRBJ1-4*01'))
    # parse_imgt_webpage(url='https://www.imgt.org/ligmdb/view.action?id=IMGT000132')