from generate_cdr3_olga import convert_olga_to_imgt_id


def read_marginals_file(filename):
    usages = {}
    read = False
    with open(filename, 'r') as infile:
        for line in infile:
            line = line.strip('\n')
            if line.startswith('@'):
                element_name = line.lstrip('@')
                if element_name in ['v_choice', 'j_choice']:
                    read = True
                else:
                    read = False
            if line.startswith('$Dim'):
                dimensions = [int(e) for e in line.lstrip('$Dim').replace('[', '').replace(']', '').split(',')]
            if line.startswith('#'):
                pass
            if line.startswith('%'):
                if read:
                    usages.update({element_name: [float(e) for e in line.lstrip('%').split(',')]})

    return usages


def get_olga_vj_usage():
    vs = {}
    js = {}
    file = '..\\venv\\Lib\\site-packages\\olga\\default_models\\mouse_T_beta\\model_marginals.txt'
    usages = read_marginals_file(file)
    v = usages['v_choice']
    j = usages['j_choice']
    for i in range(len(v)):
        try:
            imgt_id = convert_olga_to_imgt_id(i, 'V')
            vs.update({imgt_id: v[i]*100})
        except KeyError:
            print('problem, V-gene-IMGT-id not found. OLGA_id and usage:')
            print(i, v[i])

    for i in range(len(j)):
        try:
            imgt_id = convert_olga_to_imgt_id(i, 'J')
            js.update({imgt_id: j[i]*100})
        except KeyError:
            print('problem, J-gene-IMGT-id id not found. OLGA_id and usage:')
            print(i, j[i])
    print(vs)
    print(js)
    return vs, js


if __name__ == '__main__':
    get_olga_vj_usage()
