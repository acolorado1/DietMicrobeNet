import xml.etree.ElementTree as ET
import csv

parser = arg.ArgumentParser(description="Takes downloaded xml file from hdmb metabolites and subsets it to contain"
"only kegg IDs and chemical taxonomy information in CSV format")
parser.add_argument("--h_met", "-hmdb_met", type=str, help="file path for XML formatted file")
parser.add_argument("--o", "-out_file", type=str, help="files path to output CSV")
args = parser.parse_args()
def subset_hmdb_xml(hmdb_met:str, out_file:str):
    """subset of the whole hmdb metabolites file so as only to get kegg compound IDs
    and chemical taxonomy

    Args:
        hmdb_met (str): file path for downloaded hmdb xml metabolites file 
        out_file (str): output path for subsetted csv 
    """
    # Load the XML file
    tree = ET.parse(hmdb_met)
    root = tree.getroot()

    # define namespace 
    namespaces = {'ns': "http://www.hmdb.ca"}

    # Prepare the CSV file
    csv_file = out_file  # Replace with the desired output file name
    csv_fields = ["kegg_id", "taxonomy_direct_parent",
                "taxonomy_kingdom", "taxonomy_super_class", "taxonomy_class", 
                "taxonomy_sub_class"]

    # Open CSV file for writing
    with open(csv_file, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.DictWriter(file, fieldnames=csv_fields)
        writer.writeheader()
        
        # Loop through each metabolite
        for metabolite in root.findall("ns:metabolite", namespaces=namespaces):
            #kegg_id = metabolite.find("ns:kegg_id", namespaces=namespaces).text if metabolite.find("ns:kegg_id") is not None else ""
            kegg_id_elem = metabolite.find("ns:kegg_id", namespaces)
            kegg_id = kegg_id_elem.text if kegg_id_elem is not None else None

            
            # Extract taxonomy details
            taxonomy = metabolite.find("ns:taxonomy", namespaces)
            if taxonomy is not None:
                taxonomy_details = {
                    "taxonomy_direct_parent": taxonomy.find("ns:direct_parent", namespaces).text if taxonomy.find("ns:direct_parent", namespaces) is not None else "",
                    "taxonomy_kingdom": taxonomy.find("ns:kingdom", namespaces).text if taxonomy.find("ns:kingdom", namespaces) is not None else "",
                    "taxonomy_super_class": taxonomy.find("ns:super_class", namespaces).text if taxonomy.find("ns:super_class", namespaces) is not None else "",
                    "taxonomy_class": taxonomy.find("ns:class", namespaces).text if taxonomy.find("ns:class", namespaces) is not None else None,
                    "taxonomy_sub_class": taxonomy.find("ns:sub_class", namespaces).text if taxonomy.find("ns:sub_class", namespaces) is not None else "",
                }
            else:
                taxonomy_details = {key: None for key in csv_fields[1:]}

            # Combine kegg_id and taxonomy details
            row = {"kegg_id": kegg_id} | taxonomy_details
            writer.writerow(row)

    print(f"Data successfully extracted to {csv_file}.")

#subset_hmdb_xml("/Users/burkhang/Code_Projs/DietMicrobeNet/Data/hmdb_metabolites.xml", "/Data/hmdb.csv")
subset_hmdb_xml(args.h_met, args.o)