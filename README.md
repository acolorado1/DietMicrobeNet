# DietMicrobeNet: Network Modeling of Dietary Effect on Microbial Metabolism

The purpose of this code will be to create a metabolic network where nodes represent compounds and edges represent reactions. Compounds will originate from either food items (e.g. apple), microbes, both or neither. Using [Neo4j](https://neo4j.com/?utm_source=GSearch&utm_medium=PaidSearch&utm_campaign=Evergreen&utm_content=AMS-Search-SEMBrand-Evergreen-None-SEM-SEM-NonABM&utm_term=neo4j&utm_adgroup=core-brand&gad_source=1&gad_campaignid=20973570619&gbraid=0AAAAADk9OYrR1oC8MGny6LkW0dPjJNXme&gclid=CjwKCAjw6P3GBhBVEiwAJPjmLj_XgC9eELFrA3JqCbNXKo_gQfzsk-rvJUDym_R2Bo9ccBjtzoXz_BoCCLwQAvD_BwE), a graph database management system, graphs developed can be queried to find instances where microbes have the potential to metabolize dietary compounds. Information such as food frequency, read abundance, and taxonomy can also be conserved within the graph depending on user preferences.

## Install 

In the terminal, go to directory of choice and clone this repo:

```
git clone https://github.com/acolorado1/DietMicrobeNet.git       # clone repo
cd DietMicrobeNet                                                # move into this project directory 
```

Create environment with yaml file provided:

```
conda env create -f DMnet_env.yaml                              # create environment
conda activate DietMicrobeNet                                   # activate environment 
pip install -e .                                                # set up directory structure 
```

## Workflow 

### Downloading FooDB and HMDB database information

For the following scripts to run you will need four files taken from FooDB and HMDB located in a public drive. 

To do this run: 

```
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/ESXx7vpypQFOt4iVv6x-ErkBykpAVS1fppQjYZkrxkDnAA?download=1' -O Data/CompoundExternalDescriptor.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EYJUYQWmY9VDlYZIAXpzpvEBzhrnViFZQjrikXIla_aPPg?download=1' -O Data/Content.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EXyRAlYs1htNlcwz5T67BxQBGO7HfOjmfIBlkOydM0BIAw?download=1' -O Data/Food.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EbY2fD3JTcNLomKFqQhY5jABAXN-60A80PmkngRynazocg?download=1' -O Data/hmdb.csv
```

### Getting Compound Lists 

A web app has been developed to create CSVs that will reference specific food items either from KEGG or FooDB and user will input the frequency of their consumption. 

Then, for web app usage, run: 

```
# in terminal
streamlit run src/get_foods.py
```
1. Search for food 
2. Add user input value from 1 to 100 (frequency of consumption as a percentage)
3. Download filtered dataset
4. Follow shut down instructions 

### Inputs Needed for Graph Creation

To create the graph, you will need a directory containing 3-4 files. 

1. ko_taxonomy_abundance.csv 
2. noquote_ko.txt
3. foodb_foods_dataframe.csv AND/OR kegg_organisms_dataframe.csv (downloaded from streamlit app)

For example, if you are running both metabolome and whole genome methods, directory structure should look as following: 

``` bash
path/to/dir1/
├──ko_taxonomy_abundance.csv
├──noquote_ko.txt
├──foodb_foods_dataframe.csv
└──kegg_organisms_dataframe.csv
```

#### Example files 

Things to know about *ko_taxonomy_abundance.csv* 

- Should have three columns "KO", "taxonomy", and a column representing read abundance (in this case it is "Abundance_RPKs")
- ONLY the read abundance column is mutable meaning this CSV must have these EXACT column names for KO and taxonomy 
- If you do not have taxonomy or abundance information **leave the column blank**, downstream process will eliminate empty values 
```
"KO","taxonomy","Abundance_RPKs"
"K00001","g__Bifidobacterium.s__Bifidobacterium_bifidum",30.025907407
"K00001","g__Bifidobacterium.s__Bifidobacterium_longum",0
"K00001","unclassified",0
"K00002","g__Blautia.s__Blautia_obeum",41.8831170812
"K00002","g__Blautia.s__Blautia_sp",0
"K00002","g__Blautia.s__Blautia_sp_AF19_10LB",0
"K00002","unclassified",0
```

Things to know about *noquote_ko.txt*

- It is a list of KOs (KEGG gene IDs)
- There must be no quotes or commas in the file 

```
K00001
K00002
K00003
K00004
K00005
K00008
```

Both *foodb_foods_dataframe.csv* and *kegg_organisms_dataframe.csv* are the CSVs you can download through the streamlit app. Once downloaded the only thing that can be safely changed is either the removal of an entire row or adjustment of the food frequency. 

### Creating the Graph

Here is a list of required and optional arguments that can be accessed by typing **python run_workflow.py -h**:

```
usage: run_workflow.py [-h] --directories DIRECTORIES [DIRECTORIES ...] [--metabolome] [--genome] [--e-weights] [--n-weights] [--include-orgs] [--abundance-col ABUNDANCE_COL] [--cores CORES] [--profile PROFILE] [--dry-run]

Wrapper for Snakemake workflow

options:
  -h, --help            show this help message and exit
  --directories DIRECTORIES [DIRECTORIES ...]
                        List of input directories
  --metabolome          Enable metabolome analysis
  --genome              Enable genome analysis
  --e-weights           Enable E weights
  --n-weights           Enable N weights
  --include-orgs        Include organisms
  --abundance-col ABUNDANCE_COL
                        Column name for abundance
  --cores CORES         Number of cores to use
  --profile PROFILE     Snakemake profile to use
  --dry-run, -n         Perform a dry-run
```

Things to know about the arguments: 

 - directories should be separated by a space (see Example Usage #1)
 - metabolome, genome, e-weights, n-weights, and include-orgs should only be included if they are wanted (see Example Usage #1 and #2)
 - cores, profiles, and dry-run are Snakemake specific arguments (see [Snakemake's official documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html))

#### Example Usage #1

Running with two directories, metabolome, genome and only including edge weights. 

```
python run_workflow.py \                                                    
--directories /path/to/dir1 /path/to/dir2 \
--metabolome \
--genome \
--e-weights
```

#### Example Usage #2

Running with one directory, only metabolome, no node or edge weights, using multiple cores 

```
python run_workflow.py \                                                    
--directories /Users/burkhang/Desktop/Zim_snakes/F_Exp/ZIM091_A \
--metabolome \
--cores 4
```

### Outputs 

When you run the workflow you a directory will be created with the following structure: 

``` bash
my_directory/
├──ko_taxonomy_abundance.csv
├──noquote_ko.txt
├──foodb_foods_dataframe.csv
├──kegg_organisms_dataframe.csv
├──output_gen/
    ├──AMON_output/
        ├──AMON_log.txt
        ├──gene_set_1_enrichment.tsv
        ├──gene_set_2_enrichment.tsv
        ├──kegg_mapper.tsv
        ├──origin_table.tsv
        ├──enrichment_heatmap.png
        ├──venn.png
        ├──co_dict.json
        ├──ko_dict.json
        └──rn_dict.json
    ├──graph/
        ├──WG_edges_df.csv
        ├──WG_nodes_df.csv
        ├──WG_AbundanceDistribution.png
        └──WG_FoodFrequencyDistribution.png
    ├──org_KO/
        ├──a_file_for_every_food.txt
        └──joined.txt
    ├──food_item_kos.csv
    └──compound_report.html
├──output_met/
    ├──AMON_output/
        ├──AMON_log.txt
        ├──gene_set_1_enrichment.tsv
        ├──kegg_mapper.tsv
        ├──origin_table.tsv
        ├──enrichment_heatmap.png
        ├──co_dict.json
        ├──ko_dict.json
        └──rn_dict.json
    ├──graph/
        ├──M_edges_df.csv
        ├──M_nodes_df.csv
        ├──M_AbundanceDistribution.png
        └──M_FoodFrequencyDistribution.png
    ├──food_meta.csv
    └──compound_report.html
```

- output_gen and output_met are outputs created from whole genome and metabolome methods respectively
- in output_gen/org_KO file there will be as many text files as there are food items (this is represented by *a_file_for_every_food.txt* )
- graph outputs include node and edge dataframes as well as histograms of node and edge weights if those were included 
- compound_report.html gives information on the compounds found in each food 