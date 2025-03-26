# Creating A Diet-Microbe Network

The purpose of this code will be to create a custom collection of compounds representing a diet of the users choice and integrate that along with microbial genes to create a network where nodes represent compounds and edges genes.

This effort is to characterize microbial metabolism of dietary compounds

## Install 

In the terminal, go to directory of choice and clone this repo:

```
git clone https://github.com/acolorado1/DietMicrobeNet.git
```

Create environment with yaml file provided:

```
conda env create -f DMnet_env.yaml
conda activate DietMicrobeNet
```

## Workflow 

### Getting Compound Lists 

First determine whether you want to use a ***metabolomics food database*** or use ***whole genomes*** from organisms to predict compounds. 

**Note**: you can use both and a function will be provided to merge both results (TODO)

#### For metabolome: 

```
# in terminal
streamlit run FetchFood/get_FooDB_food.py 
```
1. Search for food 
2. Add *ROW NUMBER* to filtered dataset 
3. Download filtered dataset
4. Follow shut down instructions 

```
# in terminal 
Rscript CompoundQueries comp_FooDB.R \
--diet_file "path/to/file/downloaded/from/website.csv" \
--content_file "path/to/content/file/from/fooDB.csv" \
--ExDes_file "path/to/externaldescriptor/file/from/foodb.csv" \
--output_file "path/to/output/directory/output.txt"
--meta_o_file "path/to/output/directory/meta_foodb.csv"
```

#### For whole genomes: 

**Note**: This will take a large amount of time due to the need to KEGG API restrictions, if you have downloaded KEGG that will make it faster see [AMON github](https://github.com/lozuponelab/AMON)
```
# in terminal 
streamlit run FetchFood/get_KEGG_food.py
```
1. Search organism 
2. Add to dataframe 
3. Download final dataframe containing all organisms 
4. Follow shut down instructions

