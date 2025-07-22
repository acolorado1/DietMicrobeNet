# Creating A Diet-Microbe Network

The purpose of this code will be to create a custom collection of compounds representing a diet of the users choice and integrate that along with microbial genes to create a network where nodes represent compounds and edges genes.

This effort is to characterize microbial metabolism of dietary compounds.

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

A web app has been developed to create CSVs that will reference specific food items either from KEGG or FooDB and user will input the frequency of their consumption. 

```
# in terminal
streamlit run src/get_food.py
```
1. Search for food 
2. Add user input value from 1 to 100 (frequency of consumption as a percentage)
3. Download filtered dataset
4. Follow shut down instructions 
