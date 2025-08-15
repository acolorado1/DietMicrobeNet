import pandas as pd 
import json

def data_read_in(f_meta:str, 
                m_meta:str, 
                mapper:str, 
                rn_json:str, 
                abundance_column:str, 
                e_weights:bool, 
                orgs:bool):
    """read in data and convert to pandas dataframes or dictionaries 

    Args:
        f_meta (str): path to file containing food KOs and the foods their associated with, output of comp_KEGG.py
        m_meta (str): path to file containing microbial metadata (e.g., KOs, taxonomy, abundance)
        mapper (str): path to file from AMON output called kegg_mapper.tsv
        rn_json (str): path to file from AMON output called rn_dict.json
        abundance_column (str): name of the column in m_meta where abundance values are located
        e_weights (bool): wether abundance should be included 
        orgs (bool): True or false wether organisms will be included in edge information 

    Returns:
        three pandas dataframes (food meta, microbe meta, and mapper) and one dictionary (reactions)
    """
    
    # get food and microbial metadata and clean 
    food_meta = pd.read_csv(f_meta)

    microbe_meta = pd.read_csv(m_meta)
    if e_weights and orgs:
        m_meta_clean = microbe_meta.dropna(subset=['KO', 'taxonomy', abundance_column]).copy() # remove Nan
    elif orgs:
        m_meta_clean = microbe_meta.dropna(subset=['KO', 'taxonomy']).copy() # remove Nan
    elif e_weights:
        m_meta_clean = microbe_meta.dropna(subset=['KO', abundance_column]).copy() # remove Nan
    else:
        m_meta_clean = microbe_meta.dropna(subset=['KO']).copy() # remove Nan
    m_meta_clean['taxonomy'] = m_meta_clean['taxonomy'].astype(str) # convert taxonomy to string 
    print('NAs have been removed from microbe metadata')

    # get mapper and make origin column 
    # key for AMON mapper -> blue = gene set, green = both, yellow = other gene set 
    # when running AMON must make sure that microbe KOs are gene set 1
    mapper_df = pd.read_csv(mapper, sep='\t')
    mapper_df.columns = ['id', 'amon_origin']
    mapper_df['origin'] = mapper_df['amon_origin'].map({
    'blue': 'microbe',
    'green': 'both',
    'yellow': 'food'
    })

    # turn json in dictionary 
    with open(rn_json, 'r') as json_rn: 
        rxns = json.load(json_rn)

    return food_meta, m_meta_clean, mapper_df, rxns 

def make_organisms_abundance_dict(microbe_meta_clean, abundance_column:str, e_weights:bool, orgs:bool):
    """associated microbe taxonomy and abudance values with KOs 

    Args:
        microbe_meta_clean (pandas df): dataframe containing KOs, taxonomy and abundance 
        abundance_column (str): name of the column containing abundance information 
        e_weights (bool): True or false whether edge weights (abundance will be included)
        orgs (bool): True or false wether organisms will be included in edge information 
    
    Returns:
        two dictionaries: {KO:abundance} and {KO: organisms}
    """
    if e_weights:
        # create dict of KOs with abundance 
        ko_abundance_sum = microbe_meta_clean.groupby('KO', as_index=False).sum(abundance_column)
        ko_abundance_dict = ko_abundance_sum.set_index('KO')[abundance_column].to_dict()
    else:
        ko_abundance_dict = None

    if orgs:
        # create dict of KOs with list of associated organisms -> removed duplicates and sorted alphabetically
        ko_orgs_df = microbe_meta_clean.groupby('KO')['taxonomy'].agg(lambda x: ', '.join(sorted(set(x)))).reset_index()
        ko_orgs_dict = ko_orgs_df.set_index('KO')['taxonomy'].to_dict()
    else:
        ko_orgs_dict = None

    return ko_abundance_dict, ko_orgs_dict

def get_info_dicts(rxns: dict): 
    """from the reactions found through AMON create dictionaries to easily connect information 

    Args:
        rxns (dict): dictionary containing each reaction as key and a dictionary as values

    Returns:
        3 dictionaries: {compound: [KOs]}, {reaction: [KOs]} and {reaction: [[reactants], [products]]}
    """
    ko_compounds = {}     # {KO: set(compounds)}
    compound_kos = {}     # {compound: set(KOs)}
    rxn_kos = {}          # {reaction: set(KOs)}
    rxn_equation = {}     # {reaction: (reactants, products)}

    for rxn, rn_info_dict in rxns.items():
        # get equation info
        equation = rn_info_dict['EQUATION']
        rxn_equation[rxn] = equation

        # get KOs associated with reaction
        orthology = rn_info_dict['ORTHOLOGY']
        KOs = set()  
        for index in range(len(orthology)):
            ko = orthology[index][0]
            KOs.add(ko)

            # Update ko_compounds 
            ko_compounds.setdefault(ko, set()).update(equation[0] + equation[1])

        rxn_kos[rxn] = KOs 

        # For each KO in this reaction, update compound_kos for all involved compounds
        for ko in KOs:
            compounds = ko_compounds.get(ko, set())
            for compound in compounds:
                compound_kos.setdefault(compound, set()).add(ko)
        
    # Convert all sets to **sorted lists** to ensure stable ordering
    compound_kos = {k: sorted(v) for k, v in compound_kos.items()}
    rxn_kos = {k: sorted(v) for k, v in rxn_kos.items()}
    rxn_equation = {k: [list(equation[0]), list(equation[1])] for k, equation in rxn_equation.items()}

    return compound_kos, rxn_kos, rxn_equation

def get_organisms(kos:list, ko_organisms:dict): 
    """from a list of KOs get all associated organisms 

    Args:
        kos (list): list of KOs 
        ko_organisms (dict): dictionary where keys are KOs and values are organisms 

    Returns:
        list: organisms associated with KOs
    """
    # initialize list of organisms 
    organisms = []
    # for each ko get organisms
    organisms.extend(ko_organisms[ko] for ko in kos if ko in ko_organisms)
    
    # return set for no duplications 
    return list(set(organisms))

def get_abundance(kos:list, ko_abundance:dict): 
    """from a list of KOs get total summed abundance 

    Args:
        kos (list): list of KOs
        ko_abundance (dict): dictionary where keys are KOs and values are a number 

    Returns:
        float: total abundance of the list of KOs
    """
    # returning sum of abundance for all KOs
    return sum(ko_abundance[ko] for ko in kos if ko in ko_abundance)

def build_edges_df(mapper, 
                   rxn_kos:dict, 
                   rxn_equation:dict, 
                   orgs:bool, 
                   e_weights:bool, 
                   ko_organisms:dict, 
                   ko_abundance:dict):      # sourcery skip: for-append-to-extend
    """build a dataframe containing edge information

    Args:
        mapper (pandas df): dataframe from AMON mapper file
        rxn_kos (dict): keys are reactions and values a list of KOs
        rxn_equation (dict): keys are reactions and values a list of lists of reactants and products
        orgs (bool): whether organisms will be included 
        e_weights (bool): whether abundance information will be included 
        ko_organisms (dict): keys are kos and values a list of associated organisms 
        ko_abundance (dict): keys are kos and values are numbers representing their abundance 

    Returns:
        pandas df: each row represents an edge and associated information 
    """
    edges_list = []
    for rxn, kos in rxn_kos.items():
        # Filter: keep only microbial reactions
        subset = mapper[mapper['id'].isin(kos)]
        is_microbial = subset['origin'].isin(['microbe', 'both']).any()
        if not is_microbial:
            continue

        # Optional fields
        organisms = get_organisms(kos, ko_organisms) if orgs else pd.NA
        abundance = get_abundance(kos, ko_abundance) if e_weights else pd.NA

        reactants = rxn_equation[rxn][0]
        products = rxn_equation[rxn][1]

        for r in reactants:
            for p in products:
                edges_list.append({
                    'compound1': r,
                    'compound2': p,
                    'reaction': rxn,
                    'KOs': list(kos),
                    'organisms': organisms,
                    'abundance': abundance
                })

    return pd.DataFrame(edges_list)

def build_nodes_df(mapper, food_meta, compound_kos: dict, frequency: bool):
    """build a dataframe containing information about all nodes 

    Args:
        mapper (pandas df): AMON mapper file output 
        food_meta (pandas df): dataframe containing food KOs and their associated food items
        compound_kos (dict): keys are compounds and values a list of KOs 
        frequency (bool): whether food frequency will be taken into account 

    Returns:
        pandas df: each row represents a node and associated information 
    """
    nodes_list = []

    for compound in sorted(compound_kos.keys()):
        kos = compound_kos[compound]
        
        # get origin of compounds from AMON mapper file
        origin_series = mapper.loc[mapper['id'] == compound, 'origin']
        if len(origin_series) == 1:
            origin = origin_series.iloc[0]
        elif len(origin_series) == 0:
            origin = 'none'
        else:
            origin = '>1'

        # Default to NA
        foods = pd.NA
        freq = pd.NA

        # Only assign food associations if origin is food or both
        if origin in ('food', 'both'):
            subset = food_meta[food_meta['kegg_id'].isin(kos)]
            foods = list(subset['ScientificName'].unique())
            
            if frequency: # if nodes are weights by frequency of consumption 
                unique_foods = subset.drop_duplicates(subset=['ScientificName'])
                freq = unique_foods['food_frequency'].sum()

        nodes_list.append({
            'compound': compound,
            'origin': origin,
            'assoc_food': foods,
            'freq': freq
        })

    return pd.DataFrame(nodes_list)
        
def summarize_res(compound_kos:dict, rxn_kos:dict,  
                  nodes_df, edges_df): 
    """create a summary of information for the network 

    Args:
        compound_kos (dict): keys are compounds and values are KOs 
        rxn_kos (dict): keys are reactions and values ae KOs 
        nodes_df (pandas df): dataframe where each row is a node 
        edges_df (pandas df): dataframe where each row is an edge 
    """
    # how many reactions were gathered 
    print(f'{len(rxn_kos)} number of reactions were found')

    # how many KOs mapped to a reaction 
    all_unique_kos = set(sum(compound_kos.values(), []))
    print(f'{len(all_unique_kos)} KOs mapped to reactions')
    
    # how many nodes were made 
    print(f'{len(nodes_df)} nodes were created')

    # how many edges 
    print(f'{len(edges_df)} edges were created')

    # how many nodes with none, both, microbial, and food origins 
    print('Number of nodes per origin category')
    print(nodes_df['origin'].value_counts())


    




