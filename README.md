# DietMicrobeNet: Network Modeling of Dietary Effect on Microbial Metabolism

The purpose of this code will be to create a metabolic network where nodes represent compounds and edges represent reactions. Compounds will originate from either food items (e.g. apple), microbes, both or neither. Using [Neo4j](https://neo4j.com/?utm_source=GSearch&utm_medium=PaidSearch&utm_campaign=Evergreen&utm_content=AMS-Search-SEMBrand-Evergreen-None-SEM-SEM-NonABM&utm_term=neo4j&utm_adgroup=core-brand&gad_source=1&gad_campaignid=20973570619&gbraid=0AAAAADk9OYrR1oC8MGny6LkW0dPjJNXme&gclid=CjwKCAjw6P3GBhBVEiwAJPjmLj_XgC9eELFrA3JqCbNXKo_gQfzsk-rvJUDym_R2Bo9ccBjtzoXz_BoCCLwQAvD_BwE), a graph database management system, graphs developed can be queried to find instances where microbes have the potential to metabolize dietary compounds. Information such as food frequency, read abundance, and taxonomy can also be conserved within the graph depending on user preferences.

## General Workflow 

Every user will need to set up and install this program the same way. To do this go to [1. Getting Started](https://github.com/acolorado1/DietMicrobeNet/wiki/1.-Getting-Started) in the repo's [Wiki]((https://github.com/acolorado1/DietMicrobeNet/wiki)). 

After getting set up, you will have decide if you want to run each step manually or use the provided Snakemake Workflow. 

* For manual approach follow steps 2-5 in the [Wiki](https://github.com/acolorado1/DietMicrobeNet/wiki).
* For Snakemake Workflow follow [step 2](https://github.com/acolorado1/DietMicrobeNet/wiki/2.-Find-Food-Items), then run [Snakemake](https://github.com/acolorado1/DietMicrobeNet/wiki/Snakemake-Workflow-(Steps-3%E2%80%904)), then run [step 5](https://github.com/acolorado1/DietMicrobeNet/wiki/5.-Inter%E2%80%90Sample-Comparison) all found in the [Wiki](https://github.com/acolorado1/DietMicrobeNet/wiki).

## Contact

For all programmatic questions contact: 

- Angela Sofia Burkhart Colorado (main developer) at angelasofia.burkhartcolorado@cuanschutz.edu

This work was done under the auspices of the [Lozupone Lab](https://www.lozuponelab.com/). If you have any questions regarding future publications and general management of this tool contact:

- Dr. Catherine Lozupone (PI) at catherine.lozupone@cuanschutz.edu