# Creating A Diet-Microbe Network

The purpose of this code will be to create a custom collection of compounds representing a diet of the users choice and integrate that along with microbial genes to create a network where nodes represent compounds and edges genes.

This effort is to characterize microbial metabolism of dietary compounds

## Install 

In the terminal, go to directory of choice and clone this repo.

## Workflow 

### Getting Compound Lists 

First determine whether you want to use a ***metabolomics food database*** or use ***whole genomes*** from organisms to predict compounds. 

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
--output_file "path/to/output/directory.txt"
```

#### For whole genomes: 

```
# in terminal 
streamlit run FetchFood/get_KEGG_food.py
```
1. Search organism 
2. Add to dataframe 
3. Download final dataframe containing all organisms 
4. Follow shut down instructions

