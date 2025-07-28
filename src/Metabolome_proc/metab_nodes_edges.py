import json
import pandas as pd

def node_df(microbome_compounds:str, food_compounds:str, n_weights:bool):
    """create node dataframe used as input for graph visualization 

    Args:
        microbome_compounds (str): path to compound json file from AMON output 
        food_compounds (str): path to CSV output from comp_FoodDB.py 
        n_weights (bool): true or false if food frequency is to be used 
    
    Returns: 
        pandas dataframe: dataframe containing compound information including compound
        origin, id, name, foods associated, and estimated quantity as a ratio
    """
    # read in microbe json file 
    with open(microbome_compounds, 'r') as json_comp: 
        m_comp = json.load(json_comp)

    # get a list of compounds and remove glycans 
    m_comps = list(m_comp.keys())
    for m_comp in m_comps:
        if m_comp[0] == 'G':
            m_comps.remove(m_comp)
    print(f"There are {len(m_comps)} compounds assocaited with microbes")

    # get food assocaited compounds 
    f_comp_df = pd.read_csv(food_compounds)

    # get list of unique compounds from food 
    f_comps = list(set(f_comp_df['kegg_id'].tolist()))
    print(f"There are {len(f_comps)} compounds assocaited with food items")

    # get list of all unique compounds 
    all_comps = list(set(m_comps + f_comps))
    print(f"There are {len(all_comps)} unique compounds")

    # get the origin of each compound (microbe, food, or both)
    origin_dict = {}
    for compound in all_comps: 
        if compound in m_comps and compound in f_comps:
            origin_dict[compound] = 'both'
        elif compound in m_comps: 
            origin_dict[compound] = 'microbe'
        elif compound in f_comps: 
            origin_dict[compound] = 'food'

    food_item_dict = {}
    food_freq_dict = {}

    # get all foods associated with compound and frequency of consumption
    for compound in all_comps: 

        # if associated with a food
        if compound in f_comps: 

            # subset food dataframe by compound 
            comp_df = f_comp_df[f_comp_df['kegg_id'] == compound]

            # get a list of unique food items 
            foods = list(set(comp_df['name']))
            food_item_dict[compound] = foods

            # sum frequency of consumption and make sure there is no value over 100 
            if n_weights: # if weights are being used 
                freq_sub = comp_df.drop_duplicates(subset=['name'])
                freq = freq_sub['food_frequency'].sum()
                if freq > 100:
                    raise ValueError(f"The frequency of consumption is too high: {freq}")
                food_freq_dict[compound] = int(freq)
            else:
                food_freq_dict[compound] = pd.NA # default no weights is NA 
            
            

    # Convert each dict to a DataFrame
    df_origin = pd.DataFrame(origin_dict.items(), columns=['c_id', 'origin'])
    df_food = pd.DataFrame(food_item_dict.items(), columns=['c_id', 'assoc_food'])
    df_freq = pd.DataFrame(food_freq_dict.items(), columns=['c_id', 'freq'])

    return df_origin.merge(df_food, on='c_id').merge(df_freq, on='c_id')

def edge_df(rn_json:str, ko_meta:str, nodes_df, e_weights:bool):
    """creates a pandas dataframe containing edge information for graph visualization

    Args:
        rn_json (str): file path to AMON rn_dict.json file 
        ko_meta (str): file path to csv containing KO, taxonomy, and abundance information
        nodes_df (pandas dataframe): output of nodes_df function  
        e_weights (bool): true or false if edge weights are to be used 

    Returns: 
        pandas dataframe: dataframe containing KO information including reaction, associated 
        taxonomy, and abundance RPKs  
    """
    # read in rn json 
    with open(rn_json, 'r') as rjson: 
        rxns = json.load(rjson)

    # read in ko metadata and clean data if necessary 
    meta = pd.read_csv(ko_meta)
    meta_clean = meta.dropna(subset=['KO', 'taxonomy']).copy() # remove Nan
    meta_clean['taxonomy'] = meta_clean['taxonomy'].astype(str) # convert taxonomy to string 
    print('NAs have been removed from KO metadata')

    # create dict of KOs with abundance 
    ko_abundance_sum = meta_clean.groupby('KO', as_index=False).sum('Abundance_RPKs')
    ko_abundance_dict = ko_abundance_sum.set_index('KO')['Abundance_RPKs'].to_dict()

    # create dict of KOs with list of associated organisms -> removed duplicates and sorted alphabetically
    ko_orgs_df = meta_clean.groupby('KO')['taxonomy'].agg(lambda x: ', '.join(sorted(set(x)))).reset_index()
    ko_orgs_dict = ko_orgs_df.set_index('KO')['taxonomy'].to_dict()

    # get list of nodes 
    compounds = nodes_df['c_id'].tolist()

    # dictionary to be converted into pandas dataframe 
    edges_dict = {'compound1': [], 
                  'compound2': [],
                  'reaction': [],
                  'KOs': [], 
                  'organisms': [], 
                  'abundance_RPKs': []}

    # each key is a reaction id from KEGG
    for rxn in rxns.keys():
        rn_info_dict = rxns[rxn]
        equation = rn_info_dict['EQUATION']

        # get compound from reactant side of equation
        for reactant in equation[0]: 

            # get compound from product side of equation
            for product in equation[1]: 

                # ensure both product and reactant are nodes
                if reactant in compounds and product in compounds: 
                    edges_dict['compound1'].append(reactant)
                    edges_dict['compound2'].append(product)
                    edges_dict['reaction'].append(rxn)

                    # get KOs 
                    orthology = rn_info_dict['ORTHOLOGY']
                    KOs = []
                    KOs.extend(orthology[index][0] for index in range(0,len(orthology)))

                    # for each ko get all organisms and total abundance 
                    organisms = []
                    abundance = 0 
                    for KO in KOs: 
                        if KO in ko_orgs_dict.keys():
                            organisms.append(ko_orgs_dict[KO])
                        if KO in ko_abundance_dict.keys() and e_weights:
                            abundance += ko_abundance_dict[KO]
                    
                    # if edge weights are not selected for make abundance NA
                    if not e_weights:
                        abundance = pd.NA

                    # make sure there are no duplicate organisms 
                    organisms = list(set(organisms))

                    # add to edges dictionary 
                    edges_dict['KOs'].append(','.join(KOs))
                    edges_dict['organisms'].append(organisms)
                    edges_dict['abundance_RPKs'].append(abundance)

    # return pandas dataframe of edge info 
    return pd.DataFrame(data = edges_dict) 