# Creating A Diet-Microbe Network

The purpose of this code will be to create a custom collection of compounds representing a diet of the users choice and integrate that along with microbial genes to create a network where nodes represent compounds and edges genes.

This effort is to characterize microbial metabolism of dietary compounds.

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

To do this click links, download, and move to proper directory.

1. [Compound External Descriptors](https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/ESXx7vpypQFOt4iVv6x-ErkBykpAVS1fppQjYZkrxkDnAA?e=e3mdO8)
2. [Content](https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EYJUYQWmY9VDlYZIAXpzpvEBzhrnViFZQjrikXIla_aPPg?e=dhBvmL)
3. [Food](https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EXyRAlYs1htNlcwz5T67BxQBGO7HfOjmfIBlkOydM0BIAw?e=wv3EZ9)
4. [HMDB](https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EbY2fD3JTcNLomKFqQhY5jABAXN-60A80PmkngRynazocg?e=MjDghX)

To move files into proper directory, in terminal use this command: 

```
# where /Data/ is destination file 

mv path/to/downloaded/file.csv /Data/ 
```

### Getting Compound Lists 

A web app has been developed to create CSVs that will reference specific food items either from KEGG or FooDB and user will input the frequency of their consumption. 

Then, for web app usage run: 

```
# in terminal
streamlit run src/get_foods.py
```
1. Search for food 
2. Add user input value from 1 to 100 (frequency of consumption as a percentage)
3. Download filtered dataset
4. Follow shut down instructions 
