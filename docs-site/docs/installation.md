# Installation

Two general things are need for program usage: 

1. Environment setup 
2. Input file preparation 
   
## Setting up your environment 

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

### Download FooDB and HMDB database information

For the following scripts to run you will need four files taken from FooDB and HMDB located in a public drive, and one file for all food item compounds found in FooDB. 

To do this run: 

```
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/ESXx7vpypQFOt4iVv6x-ErkBykpAVS1fppQjYZkrxkDnAA?download=1' -O Data/CompoundExternalDescriptor.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EYJUYQWmY9VDlYZIAXpzpvEBzhrnViFZQjrikXIla_aPPg?download=1' -O Data/Content.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EXyRAlYs1htNlcwz5T67BxQBGO7HfOjmfIBlkOydM0BIAw?download=1' -O Data/Food.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EbY2fD3JTcNLomKFqQhY5jABAXN-60A80PmkngRynazocg?download=1' -O Data/hmdb.csv
wget 'https://olucdenver-my.sharepoint.com/:x:/g/personal/angelasofia_burkhartcolorado_cuanschutz_edu/EZ1pyHd616RFkR9zG6kenuoBhZDroHYTbaGmEfwpxFOHLg?download=1' -O Data/AllFood/food_meta.csv
```

## Needed Files for Input 

There are three types of files that are needed to run the program: 

1. File containing a list of KOs that were found in a sample 
2. File containing KO metadata
3. File containing a list of food items that combined represent diet 

### List of KOs 

* Needs to be named _noquote_ko.txt_
* File needs to have no quotes or commas 

Example of the file: 

```
K00001
K00002
K00003
K00004
K00005
K00008
```

### KO Metadata 

* Needs to be named _ko_taxonomy_abundance.csv_
* Should have three columns "KO", "taxonomy", and a column representing read abundance (in this case it is "Abundance_RPKs")
* ONLY the read abundance column is mutable meaning this CSV must have these EXACT column names for KO and taxonomy 
* If you do not have taxonomy or abundance information **leave the column blank**, downstream process will eliminate empty values 

Example of the file: 

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

### Food items 

* These are typically downloaded in the second step of this workflow (see [2. Find Food Items](https://github.com/acolorado1/DietMicrobeNet/wiki/2.-Find-Food-Items)) 
* Named either *foodb_foods_dataframe.csv* or *kegg_organisms_dataframe.csv*
* If you want to include all possible food items you can use *Data/AllFood/food_meta.csv* (see [2. Find Food Items](https://github.com/acolorado1/DietMicrobeNet/wiki/2.-Find-Food-Items)) 