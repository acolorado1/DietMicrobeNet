import os

# -----------------------------------
# Configuration
# -----------------------------------
DIRECTORIES = config["directories"].split(",")
FOODB    = config.get("foodb", False)
GENOME        = config.get("genome", False)
E_WEIGHTS     = config.get("e_weights", False)
N_WEIGHTS     = config.get("n_weights", False)
INCLUDE_ORGS  = config.get("include_orgs", False)
ABUNDANCE_COL = config.get("abundance_col", " ") 
ALL_FOOD      = config.get("all_food", False)

print("Running with config:")
print(f"  Directories:    {DIRECTORIES}")
print(f"  FooDB:          {FOODB}")
print(f"  Genome:         {GENOME}")
print(f"  E Weights:      {E_WEIGHTS}")
print(f"  N Weights:      {N_WEIGHTS}")
print(f"  Include Orgs:   {INCLUDE_ORGS}")
print(f"  Abundance Col:  {ABUNDANCE_COL}")
print(f"  All Food:       {ALL_FOOD}")

# -------------------------------------------------------------
# Select correct food_meta file depending on ALL_FOOD
# -------------------------------------------------------------
def select_meta_file(wildcards):
    """
    If ALL_FOOD == True, use precomputed metadata from Data/AllFood.
    Otherwise use standard output_fdb/food_meta.csv and require the
    CreateFoodMetadata_fdb rule to generate it.
    """
    if ALL_FOOD:
        return f"Data/AllFood/food_meta.csv"
    else:
        return f"{wildcards.dir}/output_fdb/food_meta.csv"


# -----------------------------------
# Rule all – gather outputs across dirs
# -----------------------------------
rule all:
    input:
        # Foodb requirements 
        (expand("{dir}/output_fdb/food_compound_report.html", dir=DIRECTORIES) if FOODB and not ALL_FOOD else []),
        (expand("{dir}/output_fdb/microbe_compound_report.html", dir=DIRECTORIES) if FOODB and INCLUDE_ORGS and N_WEIGHTS else []),
        (expand("{dir}/output_fdb/graph/graph_results_report.html", dir=DIRECTORIES) if FOODB else []),

        # Genome requirements
        (expand("{dir}/output_gen/food_compound_report.html", dir=DIRECTORIES) if GENOME else []),
        (expand("{dir}/output_gen/microbe_compound_report.html", dir=DIRECTORIES) if GENOME and INCLUDE_ORGS and N_WEIGHTS else []),
        (expand("{dir}/output_gen/graph/graph_results_report.html", dir=DIRECTORIES) if GENOME else [])


# ---------------------------
# FOODB rules
# ---------------------------
if FOODB:

    rule all_fdb:
        input: 
            "{dir}/output_fdb/food_compound_report.html",
            "{dir}/output_fdb/AMON_output/rn_dict.json",
            "{dir}/output_fdb/graph/M_nodes_df.csv",
            "{dir}/output_fdb/graph/M_edges_df.csv",
            "{dir}/output_fdb/graph/M_AbundanceDistribution.png",
            "{dir}/output_fdb/graph/M_FoodFrequencyDistribution.png", 
            "{dir}/output_fdb/graph/network_summary.txt",
            "{dir}/output_fdb/microbe_compound_report.html" if INCLUDE_ORGS and N_WEIGHTS else [],
            "{dir}/output_fdb/graph/graph_results.csv",
            "{dir}/output_fdb/graph/graph_results_report.html"

    rule CreateFoodMetadata_fdb:
        input: 
            f_file = "{dir}/foodb_foods_dataframe.csv"
        output: 
            f_meta = "{dir}/output_fdb/food_meta.csv"
        conda: "DMnet_env.yaml"
        shell:
            """
            Rscript src/Metabolome_proc/comp_FoodDB.R \
                --diet_file {input.f_file} \
                --content Data/Content.csv \
                --ExDes_file Data/CompoundExternalDescriptor.csv \
                --meta_o_file {output.f_meta}
            """
    if not ALL_FOOD:
        rule CreateCompoundReport_fdb:
            input: 
                f_meta = select_meta_file,
                graphs = "{dir}/output_fdb/graph/M_nodes_df.csv" 
            output: 
                report = "{dir}/output_fdb/food_compound_report.html"
            conda: "DMnet_env.yaml"
            shell:
                """
                python {workflow.basedir}/src/Metabolome_proc/RenderCompoundAnalysis.py \
                    --food_file {input.f_meta} \
                    --output {output.report}
                """
    
    rule PrepareAMONOutput_fdb:
        input: 
            dir="{dir}"
        output:
            touch("{dir}/output_fdb/AMON_output/.prepared")
        run:
            import os, shutil
            outdir = os.path.join(input.dir, "output_fdb", "AMON_output")
            if os.path.exists(outdir):
                shutil.rmtree(outdir)
            os.makedirs(outdir, exist_ok=True)
            # create a dummy file so Snakemake sees this as complete
            open(output[0], 'w').close()


    rule RunAMON_fdb:
        input: 
            prep="{dir}/output_fdb/AMON_output/.prepared",
            kos = "{dir}/noquote_ko.txt"
        output: 
            rn_json = "{dir}/output_fdb/AMON_output/rn_dict.json"
        conda: "DMnet_env.yaml"
        shell:
            """
            rm -rf {wildcards.dir}/output_fdb/AMON_output
            amon \
                -i {input.kos} \
                -o {wildcards.dir}/output_fdb/AMON_output \
                --save_entries
            """

    rule GraphCreation_fdb:
        input: 
            f_meta = select_meta_file,
            rn_json = "{dir}/output_fdb/AMON_output/rn_dict.json",
            m_meta = "{dir}/ko_taxonomy_abundance.csv"
        params:
            flags = " ".join(filter(None, [
                "--n_weights" if N_WEIGHTS else "",
                "--e_weights" if E_WEIGHTS else "",
                "--org" if INCLUDE_ORGS else ""
            ])),
            abundance = ABUNDANCE_COL,
            graph_dir = "{dir}/output_fdb/graph/"
        output:
            nodes = "{dir}/output_fdb/graph/M_nodes_df.csv",
            edges = "{dir}/output_fdb/graph/M_edges_df.csv",
            summary = "{dir}/output_fdb/graph/network_summary.txt"
        conda: "DMnet_env.yaml"
        shell:
            """
            mkdir -p {params.graph_dir}
            python src/Metabolome_proc/main_metab.py \
                --f {input.f_meta} \
                --r {input.rn_json} \
                --m_meta {input.m_meta} \
                {params.flags} \
                --a {params.abundance} \
                --o {params.graph_dir}
            """

    if INCLUDE_ORGS and N_WEIGHTS: 
        rule MicrobeCompoundReport_fdb:
            input:
                nodes = "{dir}/output_fdb/graph/M_nodes_df.csv",
                edges = "{dir}/output_fdb/graph/M_edges_df.csv"
            output:
                report = "{dir}/output_fdb/microbe_compound_report.html"
            conda: "DMnet_env.yaml"
            shell:
                """
                python {workflow.basedir}/src/RenderCompoundAnalysis_Microbe.py \
                    --node_file {input.nodes} \
                    --edge_file {input.edges} \
                    --output {output.report}
                """
    
    rule RunGraph_fdb: 
        input: 
            nodes = "{dir}/output_fdb/graph/M_nodes_df.csv",
            edges = "{dir}/output_fdb/graph/M_edges_df.csv"
        output: 
            output = "{dir}/output_fdb/graph/graph_results.csv"
        conda: "DMnet_env.yaml"
        shell:
            """
            python src/run_graph.py \
                --n {input.nodes} \
                --e {input.edges} \
                --o {output.output}
            """
        
    rule PatternReport_fdb: 
        input: 
            graph_res = "{dir}/output_fdb/graph/graph_results.csv",
            rxn_json = "{dir}/output_fdb/AMON_output/rn_dict.json"
        output: 
            output = "{dir}/output_fdb/graph/graph_results_report.html"
        conda: "DMnet_env.yaml"
        shell: 
            """
            python {workflow.basedir}/src/RenderGraphResults_Report.py \
                --patterns {input.graph_res} \
                --rxn_json {input.rxn_json} \
                --output {output.output}
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
            "{dir}/output_gen/graph/network_summary.txt",
            "{dir}/output_gen/microbe_compound_report.html" if INCLUDE_ORGS and N_WEIGHTS else [],
            "{dir}/output_gen/graph/graph_results.csv", 
            "{dir}/output_gen/graph/graph_results_report.html"

    rule CreateFoodMetadata_gen:
        input: 
            kegg_orgs = "{dir}/kegg_organisms_dataframe.csv"
        params: 
            kos_dir = "{dir}/output_gen/org_KO"
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
            amon \
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
            flags = " ".join(filter(None, [
                "--n_weights" if N_WEIGHTS else "",
                "--e_weights" if E_WEIGHTS else "",
                "--org" if INCLUDE_ORGS else ""
            ])),
            abundance = ABUNDANCE_COL,
            graph_dir = "{dir}/output_gen/graph"
        output:
            nodes = "{dir}/output_gen/graph/WG_nodes_df.csv",
            edges = "{dir}/output_gen/graph/WG_edges_df.csv",
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
                {params.flags} \
                --a {params.abundance} \
                --o {params.graph_dir}
            """

    
    rule CreateCompoundReport_gen:
        input: 
            node_file = "{dir}/output_gen/graph/WG_nodes_df.csv"
        output: 
            o = "{dir}/output_gen/food_compound_report.html"
        shell: 
            """
            python {workflow.basedir}/src/WholeGenome_proc/RenderCompoundAnalysis.py \
                --node_file {input.node_file} \
                --output {output.o}
            """
    
    if INCLUDE_ORGS and N_WEIGHTS:
        rule MicrobeCompoundReport_gen:
            input:
                nodes = "{dir}/output_gen/graph/WG_nodes_df.csv",
                edges = "{dir}/output_gen/graph/WG_edges_df.csv"
            output:
                report = "{dir}/output_gen/microbe_compound_report.html"
            conda: "DMnet_env.yaml"
            shell:
                """
                python {workflow.basedir}/src/RenderCompoundAnalysis_Microbe.py \
                    --node_file {input.nodes} \
                    --edge_file {input.edges} \
                    --output {output.report}
                """
    
    rule RunGraph_gen: 
        input: 
            nodes = "{dir}/output_gen/graph/WG_nodes_df.csv",
            edges = "{dir}/output_gen/graph/WG_edges_df.csv"
        output: 
            output = "{dir}/output_gen/graph/graph_results.csv"
        conda: "DMnet_env.yaml"
        shell:
            """
            python src/run_graph.py \
                --n {input.nodes} \
                --e {input.edges} \
                --o {output.output}
            """
    
    rule PatternReport_gen: 
        input: 
            graph_res = "{dir}/output_gen/graph/graph_results.csv",
            rxn_json = "{dir}/output_gen/AMON_output/rn_dict.json"
        output: 
            output = "{dir}/output_gen/graph/graph_results_report.html"
        conda: "DMnet_env.yaml"
        shell: 
            """
            python {workflow.basedir}/src/RenderGraphResults_Report.py \
                --patterns {input.graph_res} \
                --rxn_json {input.rxn_json} \
                --output {output.output}
            """