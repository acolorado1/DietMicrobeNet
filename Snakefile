import os

# -----------------------------------
# Configuration
# -----------------------------------
DIRECTORIES   = config.get("directories", [])
METABOLOME    = config.get("metabolome", False)
GENOME        = config.get("genome", False)
E_WEIGHTS     = config.get("e_weights", False)
N_WEIGHTS     = config.get("n_weights", False)
INCLUDE_ORGS  = config.get("include_orgs", False)
ABUNDANCE_COL = config.get("abundance_col", " ")

print("Running with config:")
print(f"  Directories:   {DIRECTORIES}")
print(f"  Metabolome:    {METABOLOME}")
print(f"  Genome:        {GENOME}")
print(f"  E Weights:     {E_WEIGHTS}")
print(f"  N Weights:     {N_WEIGHTS}")
print(f"  Include Orgs:  {INCLUDE_ORGS}")
print(f"  Abundance Col: {ABUNDANCE_COL}")

# -----------------------------------
# Rule all – gather outputs across dirs
# -----------------------------------
rule all:
    input:
        (expand("{dir}/output_met/graph/network_summary.txt", dir=DIRECTORIES) if METABOLOME else []),
        (expand("{dir}/output_gen/graph/network_summary.txt", dir=DIRECTORIES) if GENOME else [])

# ---------------------------
# Metabolome rules
# ---------------------------
if METABOLOME:

    rule all_met:
        input: 
            "{dir}/output_met/food_meta.csv",
            "{dir}/output_met/food_compound_report.html",
            "{dir}/output_met/AMON_output/rn_dict.json",
            "{dir}/output_met/graph/M_nodes_df.csv",
            "{dir}/output_met/graph/M_edges_df.csv",
            "{dir}/output_met/graph/M_AbundanceDistribution.png",
            "{dir}/output_met/graph/M_FoodFrequencyDistribution.png", 
            "{dir}/output_met/graph/network_summary.txt"

    rule CreateFoodMetadata_met:
        input: f_file = "{dir}/foodb_foods_dataframe.csv"
        output: f_meta = "{dir}/output_met/food_meta.csv"
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
        input: 
            f_meta = "{dir}/output_met/food_meta.csv",
            graphs = "{dir}/output_met/graph/M_nodes_df.csv" 
        output: report = "{dir}/output_met/food_compound_report.html"
        conda: "DMnet_env.yaml"
        shell:
            """
            python src/Metabolome_proc/RenderCompoundAnalysis.py \
                --food_file {input.f_meta} \
                --output {output.report}
            """
    
    rule PrepareAMONOutput_met:
        input: 
            dir="{dir}"
        output:
            touch("{dir}/output_met/AMON_output/.prepared")
        run:
            import os, shutil
            outdir = os.path.join(input.dir, "output_met", "AMON_output")
            if os.path.exists(outdir):
                shutil.rmtree(outdir)
            os.makedirs(outdir, exist_ok=True)
            # create a dummy file so Snakemake sees this as complete
            open(output[0], 'w').close()


    rule RunAMON_met:
        input: 
            prep="{dir}/output_met/AMON_output/.prepared",
            kos = "{dir}/noquote_ko.txt"
        output: rn_json = "{dir}/output_met/AMON_output/rn_dict.json"
        conda: "DMnet_env.yaml"
        shell:
            """
            rm -rf {wildcards.dir}/output_met/AMON_output
            amon.py \
                -i {input.kos} \
                -o {wildcards.dir}/output_met/AMON_output \
                --save_entries
            """


    rule GraphCreation_met:
        input: 
            f_meta = "{dir}/output_met/food_meta.csv",
            rn_json = "{dir}/output_met/AMON_output/rn_dict.json",
            m_meta = "{dir}/ko_taxonomy_abundance.csv"
        params:
            n_weights_flag = "--n_weights" if N_WEIGHTS else "",
            e_weights_flag = "--e_weights" if E_WEIGHTS else "",
            org_flag = "--org" if INCLUDE_ORGS else "",
            abundance = ABUNDANCE_COL,
            graph_dir = "{dir}/output_met/graph/"
        output:
            nodes = "{dir}/output_met/graph/M_nodes_df.csv",
            edges = "{dir}/output_met/graph/M_edges_df.csv",
            a_dis = "{dir}/output_met/graph/M_AbundanceDistribution.png",
            f_dis = "{dir}/output_met/graph/M_FoodFrequencyDistribution.png",
            summary = "{dir}/output_met/graph/network_summary.txt"
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
# Genome rules
# ---------------------------
if GENOME:

    rule all_gen:
        input: 
            "{dir}/output_gen/food_item_kos.csv",
            "{dir}/output_gen/org_KO/joined.txt",
            "{dir}/output_gen/AMON_output/rn_dict.json",
            "{dir}/output_gen/AMON_output/kegg_mapper.tsv", 
            "{dir}/output_gen/graph/WG_nodes_df.csv",
            "{dir}/output_gen/graph/WG_edges_df.csv",
            "{dir}/output_gen/graph/WG_AbundanceDistribution.png",
            "{dir}/output_gen/graph/WG_FoodFrequencyDistribution.png", 
            "{dir}/output_gen/food_compound_report.html",
            "{dir}/output_gen/graph/network_summary.txt"

    rule CreateFoodMetadata_gen:
        input: kegg_orgs = "{dir}/kegg_organisms_dataframe.csv"
        params: kos_dir = "{dir}/output_gen/org_KO/"
        output: 
            food_meta = "{dir}/output_gen/food_item_kos.csv",
            joined = "{dir}/output_gen/org_KO/joined.txt"
        conda: "DMnet_env.yaml"
        shell:
            """
            mkdir -p {params.kos_dir}
            python src/WholeGenome_proc/comp_KEGG.py \
                -i {input.kegg_orgs} \
                -k {params.kos_dir} \
                -o {output.food_meta}
            """

    rule PrepareAMONOutput_gen:
        input: 
            dir="{dir}"
        output:
            touch("{dir}/output_gen/AMON_output/.prepared")
        run:
            import os, shutil
            outdir = os.path.join(input.dir, "output_gen", "AMON_output")
            if os.path.exists(outdir):
                shutil.rmtree(outdir)
            os.makedirs(outdir, exist_ok=True)
            # create a dummy file so Snakemake sees this as complete
            open(output[0], 'w').close()

    rule RunAMON_gen:
        input: 
            prep = "{dir}/output_gen/AMON_output/.prepared",  # depends on folder prep
            microbe_kos = "{dir}/noquote_ko.txt",
            diet_kos    = "{dir}/output_gen/org_KO/joined.txt"
        output: 
            rn_json = "{dir}/output_gen/AMON_output/rn_dict.json",
            mapper  = "{dir}/output_gen/AMON_output/kegg_mapper.tsv"
        conda: "DMnet_env.yaml"
        shell:
            """
            rm -rf {wildcards.dir}/output_gen/AMON_output
            amon.py \
                -i {input.microbe_kos} \
                -o {wildcards.dir}/output_gen/AMON_output \
                --other_gene_set {input.diet_kos} \
                --save_entries
            """
    
    rule GraphCreation_gen: 
        input: 
            f_meta = "{dir}/output_gen/food_item_kos.csv",
            m_meta = "{dir}/ko_taxonomy_abundance.csv",
            mapper = "{dir}/output_gen/AMON_output/kegg_mapper.tsv", 
            rn_json = "{dir}/output_gen/AMON_output/rn_dict.json"
        params: 
            n_weights_flag = "--n_weights" if N_WEIGHTS else "",
            e_weights_flag = "--e_weights" if E_WEIGHTS else "",
            org_flag = "--org" if INCLUDE_ORGS else "",
            abundance = ABUNDANCE_COL,
            graph_dir = "{dir}/output_gen/graph/"
        output:
            nodes = "{dir}/output_gen/graph/WG_nodes_df.csv",
            edges = "{dir}/output_gen/graph/WG_edges_df.csv",
            a_dis = "{dir}/output_gen/graph/WG_AbundanceDistribution.png",
            f_dis = "{dir}/output_gen/graph/WG_FoodFrequencyDistribution.png",
            summary = "{dir}/output_gen/graph/network_summary.txt"
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
        input: node_file = "{dir}/output_gen/graph/WG_nodes_df.csv"
        output: o = "{dir}/output_gen/food_compound_report.html"
        shell: 
            """
            python src/WholeGenome_proc/RenderCompoundAnalysis.py \
                --node_file {input.node_file} \
                --output {output.o}
            """