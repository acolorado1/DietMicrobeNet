# DietMicrobeNet: Network Modeling of Dietary Effect on Microbial Metabolism

The purpose of this code will be to create a metabolic network where nodes represent compounds and edges represent reactions. Compounds will originate from either food items (e.g. apple), microbes, both or neither. Graphs are created in memory, where patterns of origin (e.g., where a food originating compound is connected to a microbe originating compound) will be found and stored. Information such as food frequency, read abundance, and taxonomy can also be conserved within the graph depending on user preferences. This will allow a user to find potential instances of dietary metabolism by microbes.

## General Workflow 

Every user will need to set up and install this program the same way. To do this go to [1. Getting Started](https://github.com/acolorado1/DietMicrobeNet/wiki/1.-Getting-Started) in the repo's [Wiki](https://github.com/acolorado1/DietMicrobeNet/wiki). 

After getting set up, you will have decide if you want to run each step manually or use the provided Snakemake Workflow. 

* For **manual approach** follow steps 2-6 in the [Wiki](https://github.com/acolorado1/DietMicrobeNet/wiki).
* For **Snakemake Workflow** (RECOMMENDED) follow [step 2](https://github.com/acolorado1/DietMicrobeNet/wiki/2.-Find-Food-Items), then run [Snakemake](https://github.com/acolorado1/DietMicrobeNet/wiki/Snakemake-Workflow-(Steps-3%E2%80%904)), then run [step 5](https://github.com/acolorado1/DietMicrobeNet/wiki/5.-Inter%E2%80%90Sample-Comparison) for graph comparison, and then run [step 6](https://github.com/acolorado1/DietMicrobeNet/wiki/Metabolome-Comparison) for metabolome comparison. All steps are found in the [Wiki](https://github.com/acolorado1/DietMicrobeNet/wiki)

## Outputs

The main outputs of this program will be: 

1. HTML report on what we know about the compounds originating from the **food** (IF using specific FFQs, if all foods are used no report is created) 
2. HTML report on what we know about the compounds originating from the **microbes** (IF taxonomic and abundance information is provided)
3. HTML report on compounds and reactions involved in dietary metabolism 
4. Figures describing distribution of node and edge weights (representing food frequency and read abundance respectively)
5. Summary text file for each graph 

If you run the graph comparison script you will generate: 

1. Three summary files on the results of three different patterns that were searched for in the graph 
2. Two figures for each pattern (6 total) describing clustered similarity of graphs 

If you include a metabolome: 

1. HTML report on which compounds found in the patterns were also found in a metabolomics experiment that you have provided. 

## Contact

For all programmatic questions contact: 

- Angela Sofia Burkhart Colorado (main developer) at angelasofia.burkhartcolorado@cuanschutz.edu

This work was done under the auspices of the [Lozupone Lab](https://www.lozuponelab.com/). If you have any questions regarding future publications and general management of this tool contact:

- Dr. Catherine Lozupone (PI) at catherine.lozupone@cuanschutz.edu