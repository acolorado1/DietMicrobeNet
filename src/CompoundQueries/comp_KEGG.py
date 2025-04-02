import subprocess as sp
import pandas as pd 
import os
import glob
import json

# TODO: make parser for both functions????

def get_KEGG_compounds(org_data:str):
    """Takes in csv generated through streamlit website of organism information 
    and uses KEGG and AMON to get list of compounds 

    Args:
        org_data (str): path to CSV file with organism info
    """
    # get food organism data frame 
    orgs = pd.read_csv(org_data)

    # for each organism get all KOs and call AMON for compounds 
    org_codes = orgs['Code'].tolist()
    for code in org_codes: 

        print("Getting info on organism: " + code )
        # write where KOs are going 
        kofile = "Data/" + code + ".txt"

        # check if path exists so as not to run multiple calls 
        if os.path.exists(kofile):
            print('path exists skipping to compound prediction for this organism')
        else:
            # get KOs for organism
            command = ["extract_ko_genome_from_organism.py", 
                    "-i", code, 
                    "-o", kofile]
            sp.run(command)

        outputfile = "Data/AMON_outputs/" + code + "/"

        # check if path exists so as not to run multiple calls 
        if os.path.exists(outputfile):
            print('path exists skipping to next organism')
            continue 
        
        # call AMON for compound prediction
        command = ["AMON.py", 
                "-i", kofile, 
                "-o", outputfile, 
                "--save_entries"]
        sp.run(command)

get_KEGG_compounds("Data/test_FooDB_foods/kegg_organisms_dataframe.csv")

# TODO: function for both file outputs of compounds 
def make_compound_outputs(AMON_outputs:str): 

    # get list of organism directories
    directories = glob.glob(f'{AMON_outputs}/*')
    compounds = {}

    # TODO: make dict of list of compounds for each organism 
    for directory in directories: 
        data = json.load(directory + 'co_dict.json')
        path = directory.split("/")
        org = path[-1]
        compounds[org] = data.keys()
    
    # TODO: delete multiple copies of compounds, retain what organism it came from 

    # TODO: return compound file, and compound info file 

