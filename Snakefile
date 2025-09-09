DIRECTORY = '/Users/burkhang/Desktop/snaketest'
FOODB_FOOD_FILE = DIRECTORY + '/foodb_foods_dataframe.csv'
KEGG_FOOD_FILE = DIRECTORY + '/kegg_organisms_dataframe.csv'
ORG_KO_FILE = DIRECTORY + '/noquote_ag_sample.txt'
ORG_META_FILE = DIRECTORY + '/ko_taxonomy_abundance.csv'
METABOLOME = False
GENOME = True
E_WEIGHTS = True
N_WEIGHTS = True
INCLUDE_ORGS = True
ABUNDANCE_COL = 'Abundance_RPKs'

import os

if METABOLOME:
    met_path = DIRECTORY + '/output_met'
    os.makedirs(met_path, exist_ok=True)
if GENOME:
    gen_path = DIRECTORY + '/output_gen'
    os.makedirs(gen_path, exist_ok=True)

# ---------------------------
# Metabolome rules
# ---------------------------
if METABOLOME:
    rule all_met:
        input: 
            met_path + "/food_meta.csv",
            met_path + "/compound_report.html",
            met_path + "/AMON_output/rn_dict.json",
            met_path + "/graph/M_nodes_df.csv",
            met_path + "/graph/M_edges_df.csv",
            met_path + "/graph/M_AbundanceDistribution.png",
            met_path + "/graph/M_FoodFrequencyDistribution.png"

    rule CreateFoodMetadata_met:
        input: f_file = FOODB_FOOD_FILE
        output: f_meta = met_path + "/food_meta.csv"
        conda: "DMnet_env.yaml"
        shell:
            """
            Rscript src/Metabolome_proc/comp_FoodDB.R \
                --diet_file {input.f_file} \
                --content Data/Content.csv \
                --ExDes_file Data/CompoundExternalDescriptor.csv \
                --meta_o_file {output.f_meta}
            """

    rule CreateCompoundReport_met:
        input: f_meta = met_path + "/food_meta.csv"
        output: report = met_path + "/compound_report.html"
        conda: "DMnet_env.yaml"
        shell:
            """
            python src/Metabolome_proc/RenderCompoundAnalysis.py \
                --food_file {input.f_meta} \
                --output {output.report}
            """

    rule RunAMON_met:
        input: kos = ORG_KO_FILE
        output: rn_json = met_path + "/AMON_output/rn_dict.json"
        conda: "DMnet_env.yaml"
        shell:
            """
            amon.py \
                -i {input.kos} \
                -o {os.path.dirname(output.rn_json)} \
                --save_entries
            """

    rule GraphCreation_met:
        input: 
            f_meta = met_path + "/food_meta.csv",
            rn_json = met_path + "/AMON_output/rn_dict.json",
            m_meta = ORG_META_FILE
        params:
            n_weights_flag = "--n_weights" if N_WEIGHTS else "",
            e_weights_flag = "--e_weights" if E_WEIGHTS else "",
            org_flag = "--org" if INCLUDE_ORGS else "",
            abundance = ABUNDANCE_COL,
            graph_dir = met_path + "/graph/"
        output:
            nodes = met_path + "/graph/M_nodes_df.csv",
            edges = met_path + "/graph/M_edges_df.csv",
            a_dis = met_path + "/graph/M_AbundanceDistribution.png",
            f_dis = met_path + "/graph/M_FoodFrequencyDistribution.png"
        conda: "DMnet_env.yaml"
        shell:
            """
            mkdir -p {params.graph_dir}
            python src/Metabolome_proc/main_metab.py \
                --f {input.f_meta} \
                --r {input.rn_json} \
                --m_meta {input.m_meta} \
                {params.n_weights_flag} \
                {params.e_weights_flag} \
                {params.org_flag} \
                --a {params.abundance} \
                --o {params.graph_dir}
            """

# ---------------------------
# GENOME rules
# ---------------------------
if GENOME:

    gen_path = DIRECTORY + '/output_gen'
    os.makedirs(gen_path, exist_ok=True)

    rule all_gen:
        input: 
            gen_path + "/food_item_kos.csv",
            gen_path + "/org_KO/joined.txt",
            gen_path + "/AMON_output/rn_dict.json",
            gen_path + "/AMON_output/kegg_mapper.tsv", 
            gen_path + "/graph/WG_nodes_df.csv",
            gen_path + "/graph/WG_edges_df.csv",
            gen_path + "/graph/WG_AbundanceDistribution.png",
            gen_path + "/graph/WG_FoodFrequencyDistribution.png", 
            gen_path + "/compound_report.html"

    rule CreateFoodMetadata_gen:
        input: kegg_orgs = KEGG_FOOD_FILE
        params: kos_dir = gen_path + '/org_KO/'
        output: 
            food_meta = gen_path + '/food_item_kos.csv',
            joined = gen_path + '/org_KO/joined.txt'
        conda: "DMnet_env.yaml"
        shell:
            """
            mkdir -p {params.kos_dir}
            python src/WholeGenome_proc/comp_KEGG.py \
                -i {input.kegg_orgs} \
                -k {params.kos_dir} \
                -o {output.food_meta}
            """

    rule RunAMON_gen:
        input: 
            microbe_kos = ORG_KO_FILE,
            diet_kos = gen_path + "/org_KO/joined.txt"
        output: 
            rn_json = gen_path + "/AMON_output/rn_dict.json",
            mapper = gen_path + "/AMON_output/kegg_mapper.tsv"
        conda: "DMnet_env.yaml"
        shell:
            """
            amon.py \
                -i {input.microbe_kos} \
                -o {gen_path}/AMON_output \
                --other_gene_set {input.diet_kos} \
                --save_entries
            """
    
    rule GraphCreation_gen: 
        input: 
            f_meta = gen_path + "/food_item_kos.csv",
            m_meta = ORG_META_FILE,
            mapper = gen_path + "/AMON_output/kegg_mapper.tsv", 
            rn_json = gen_path + "/AMON_output/rn_dict.json"
        params: 
            n_weights_flag = "--n_weights" if N_WEIGHTS else "",
            e_weights_flag = "--e_weights" if E_WEIGHTS else "",
            org_flag = "--org" if INCLUDE_ORGS else "",
            abundance = ABUNDANCE_COL,
            graph_dir = gen_path + "/graph/"
        output:
            nodes = gen_path + "/graph/WG_nodes_df.csv",
            edges = gen_path + "/graph/WG_edges_df.csv",
            a_dis = gen_path + "/graph/WG_AbundanceDistribution.png",
            f_dis = gen_path + "/graph/WG_FoodFrequencyDistribution.png"
        conda: "DMnet_env.yaml"
        shell: 
            """
            mkdir -p {params.graph_dir}
            
            python src/WholeGenome_proc/main_geno.py \
                --f_meta {input.f_meta} \
                --m_meta {input.m_meta} \
                --mapper {input.mapper} \
                --rn_json {input.rn_json} \
                {params.n_weights_flag} \
                {params.e_weights_flag} \
                --a {params.abundance} \
                {params.org_flag} \
                --o {params.graph_dir}
            """
    
    rule CreateCompoundReport_gen:
        input: node_file = gen_path + "/graph/WG_nodes_df.csv",
        output: o = gen_path + "/compound_report.html"
        shell: 
            """
            python src/WholeGenome_proc/RenderCompoundAnalysis.py \
                --node_file {input.node_file} \
                --output {output.o}
            """
    
# ---------------------------
# Master all
# ---------------------------
rule all:
    input:
        (rules.all_met.input if METABOLOME else []),
        (rules.all_gen.input if GENOME else [])
