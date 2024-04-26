from olga import load_model, generation_probability as pgen, sequence_generation as seq_gen

params_filename = '../../venv/Lib/site-packages/olga/default_models/mouse_T_beta/model_params.txt'
marginals_filename = '../../venv/Lib/site-packages/olga/default_models/mouse_T_beta/model_marginals.txt'
v_anchor_pos_file = '../../venv/Lib/site-packages/olga/default_models/mouse_T_beta/V_gene_CDR3_anchors.csv'
j_anchor_pos_file = '../../venv/Lib/site-packages/olga/default_models/mouse_T_beta/J_gene_CDR3_anchors.csv'

genomic_data = load_model.GenomicDataVDJ()
genomic_data.load_igor_genomic_data(params_filename, v_anchor_pos_file, j_anchor_pos_file)

generative_model = load_model.GenerativeModelVDJ()
generative_model.load_and_process_igor_model(marginals_filename)

pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
p = pgen_model.compute_nt_CDR3_pgen('TGTGCAAGCAGCTTAGACAGTCAAAACACCTTGTACTTT')
print(p)

seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)
print('Sequences without insertions')
for i in range(10):
    sequence = seq_gen_model.gen_rnd_prod_noins_CDR3()
    print(sequence)
print('Sequences with insertions')
for i in range(10):
    sequence = seq_gen_model.gen_rnd_prod_CDR3()
    print(sequence)
