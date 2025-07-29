import shutil
import subprocess as sp
import pandas as pd 
import os
from itertools import chain
import argparse
def get_KEGG_KOs(org_data: str, output: str, overwrite=False):  
    if os.path.exists(output):
        if overwrite:
            print(f"Overwriting existing folder: {output}")
            shutil.rmtree(output)
            os.makedirs(output)
        else:
            print(f"Folder '{output}' exists. Skipping KO retrieval.")
            return
    else:
        os.makedirs(output)

    orgs = pd.read_csv(org_data)
    org_codes = orgs['Code'].tolist()

    for code in org_codes:
        print(f"Getting info on organism: {code}")
        kofile = os.path.join(output, f"{code}.txt")

        if os.path.exists(kofile):
            print(f"{kofile} exists. Skipping.")
            continue

        result = sp.run(["extract_ko_genome_from_organism.py", "-i", code, "-o", kofile],
                        capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error fetching {code}: {result.stderr}")


def merge_organism_KOs(org_dir:str, met_output:str):
    """merge organism KOs 

    Args:
        org_dir (str): directory containing all organism files 
        met_output (str): path to streamlit CSV 
    
    Returns: 
        joined text file containing all KOs for AMON 
        CSV of KOs and codes 
    """
    # get full file paths for organism KOs 
    file_paths = []
    for root, _, files in os.walk(org_dir):
        for file in files:
            full_path = os.path.join(root, file)
            file_paths.append(full_path)

    # create empty dict for org an KOs
    org_ko = {}

    # for each set of organism KOs
    for organism in file_paths:

        KOs = []

        # open file and add KOs to list 
        with open(organism, 'r') as file:
            KOs.extend([line.strip() for line in file])
            #get organism code and add to dict 
            code = os.path.splitext(os.path.basename(organism))[0]
            org_ko[code] = KOs
    
    # create pandas dataframe for compounds metadata and write csv
    meta_long = pd.DataFrame([(key, kegg_id) for key, values in org_ko.items() for kegg_id in values], 
                             columns=['code', 'kegg_id'])
    
    # save metadata file 
    meta_long.to_csv(met_output, index=False)

    # Combine all KO entries into a single list with no duplicates 
    combined_ko_list = list(set(chain.from_iterable(org_ko.values())))
    
    # Write unique KO entries to the output file
    ko_out = org_dir + "/joined.txt"
    with open(ko_out, 'w') as file:
        file.write("\n".join(combined_ko_list))

def main():
    parser = argparse.ArgumentParser(
        description="Retrieve KEGG KOs for organisms and merge into a single file for AMON. CSV containig organism codes and KOs will also be generated"
    )

    parser.add_argument(
        "-i", "--org_info",
        type=str,
        required=True,
        help="Path to CSV file with organism info (must include 'Code' column)."
    )
    parser.add_argument(
        "-k", "--ko_output",
        type=str,
        required=True,
        help="Directory where KO files for each organism will be stored."
    )
    parser.add_argument(
        "-o", "--org_ko_output",
        type=str,
        required=True,
        help="Path to CSV file where merged organism-KO mapping will be saved."
    )

    args = parser.parse_args()

    # Run main function
    get_KEGG_KOs(org_data=args.org_info, 
                 output=args.ko_output)
    
    merge_organism_KOs(org_dir=args.ko_output, 
                       met_output=args.org_ko_output)

if __name__ == "__main__":
    main()

