import subprocess as sp
import pandas as pd 
import os
def get_KEGG_KOs(org_data:str, metadata_o:str):  
    """Takes in csv generated through streamlit website of organism information 
    and uses queries KEGG for compounds associated with organism KOs and concatenates all of them 

    Args:
        org_data (str): path to CSV file with organism info

    Returns: 
        Returns two files: list of upduplicated KOs and csv of KOs and metadata (e.g., org of origin)
    """
    # get food organism data frame 
    orgs = pd.read_csv(org_data)

    # for each organism get all KOs 
    org_codes = orgs['Code'].tolist()

    # create dict of organism and associated KOs
    org_KOs = {}

    # for each organism 
    for code in org_codes: 
        print("Getting info on organism: " + code )
        
        # write where KOs are going 
        kofile = "Data/kegg_orgs/" + code + ".txt"

        # check if path exists so as not to run multiple calls 
        if os.path.exists(kofile):
            print('path exists skipping to compound prediction for this organism')
        else:
            # get KOs for organism
            command = ["extract_ko_genome_from_organism.py", 
                    "-i", code, 
                    "-o", kofile]
            sp.run(command)

        #with open(kofile, 'r'):
def run_AMON(micro_KOs: str, diet_KOs: str, output_path: str, additional_args: list = None):
    """Takes a file containing microbiome-associated KOs and diet-associated KOs and runs AMON.

    Args:
        micro_KOs (str): Path to microbial KO file.
        diet_KOs (str): Path to diet KO file.
        output_path (str): Path for AMON output.
        additional_args (list): Optional additional arguments for AMON.
    """
    # Validate input paths
    if not os.path.exists(micro_KOs):
        raise FileNotFoundError(f"Microbial KO file not found: {micro_KOs}")
    if not os.path.exists(diet_KOs):
        raise FileNotFoundError(f"Diet KO file not found: {diet_KOs}")
    
    # Construct command
    command = ["AMON.py", 
               "-i", micro_KOs, 
               "--other_gene_set", 
               diet_KOs, "-o", 
               output_path, "--save_entries"]
    
    # If there are additional arguments add 
    if additional_args:
        command.extend(additional_args)

    print(f"Running AMON with microbiome KOs: {micro_KOs} and diet KOs: {diet_KOs}")
    try:
        result = sp.run(command, check=True, capture_output=True, text=True)
        print(f"AMON completed successfully. Output saved to: {output_path}")
        return result
    except sp.CalledProcessError as e:
        print(f"Error running AMON: {e}")
        return None

#get_KEGG_compounds("/Users/burkhang/Code_Projs/DietMicrobeNet/Data/test/app_ouput/kegg_organisms_dataframe.csv")

