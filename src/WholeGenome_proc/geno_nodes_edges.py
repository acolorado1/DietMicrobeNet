import pandas as pd 
import json
def node_edge_dfs(f_meta:str, 
                  m_meta: str, 
                  mapper:str, 
                  co_json:str, 
                  rn_json:str,
                  n_weights:bool, 
                  e_weights:bool, 
                  o_nodes:str, 
                  o_edges:str):
    """create pandas dataframe of nodes originating from whole genome predictions 

    Args:
        f_meta (str): file path to CSV created from comp_KEGG.py (food info)
        m_meta (str): file path to CSV containing KO, taxonomic, and abundance info
        mapper (str): mapper file path generated when running AMON
        co_json (str): JSON file path generated on compound information from AMON 
        rn_json (str): JSON file path generated on reaction information from AMON
        n_weights (bool): True or False whether nodes should have frequencies associated with them
        e_weights (bool): True or False whether edges should have frequencies associated with them
        o_nodes (str): directory for output nodes CSV file 
        o_edges(str): directory for output edges CSV file 

    Returns:
        two CSVs containing information for each node and edge to be used in graph visualization 
    """
    # get pandas dataframe of organism codes, KOs, and food info
    food_meta = pd.read_csv(f_meta)

    # get metadata (KOs taxonomic )
    microbe_meta = pd.read_csv(m_meta)
    m_meta_clean = microbe_meta.dropna(subset=['KO', 'taxonomy']).copy() # remove Nan
    m_meta_clean['taxonomy'] = m_meta_clean['taxonomy'].astype(str) # convert taxonomy to string 
    print('NAs have been removed from microbe metadata')

    # get origin prediction from AMON and make key for graph 
    # key for AMON mapper -> blue = gene set, green = both, yellow = other gene set 
    # when running AMON must make sure that microbe KOs are gene set 1
    mapper_df = pd.read_csv(mapper, sep='\t')
    mapper_df.columns = ['id', 'amon_origin']
    mapper_df['origin'] = mapper_df['amon_origin'].map({
    'blue': 'microbe',
    'green': 'both',
    'yellow': 'food'
    })

    # read in co_dict json file from AMON results 
    with open(co_json, 'r') as json_comp: 
        compounds = json.load(json_comp)

    # read in rn_dict json file from AMON results 
    with open(rn_json, 'r') as json_rn: 
        rxns = json.load(json_rn)

    # get list of compound ids 
    compound_ids = list(compounds.keys())

    # create dict of KOs with abundance 
    ko_abundance_sum = m_meta_clean.groupby('KO', as_index=False).sum('Abundance_RPKs')
    ko_abundance_dict = ko_abundance_sum.set_index('KO')['Abundance_RPKs'].to_dict()

    # create dict of KOs with list of associated organisms -> removed duplicates and sorted alphabetically
    ko_orgs_df = m_meta_clean.groupby('KO')['taxonomy'].agg(lambda x: ', '.join(sorted(set(x)))).reset_index()
    ko_orgs_dict = ko_orgs_df.set_index('KO')['taxonomy'].to_dict()

    # dictionary of edges information  
    edges_dict = {'compound1': [], 
                  'compound2': [],
                  'reaction': [],
                  'KOs': [], 
                  'organisms': [], 
                  'abundance_RPKs': []}

    # get set of microbe KOs 
    m_set = set(m_meta_clean['KO'].to_list())

    # compounds and food KOs 
    comp_food = {}

    # each key is a reaction id from KEGG
    for rxn in rxns.keys():
        # get all info on reaction 
        rn_info_dict = rxns[rxn]

        # get KOs associated with reaction 
        orthology = rn_info_dict['ORTHOLOGY']
        KOs = []
        KOs.extend(orthology[index][0] for index in range(len(orthology)))
        k_set = set(KOs)

        # get equation associated with reaction 
        equation = rn_info_dict['EQUATION']

        # make sure reaction is microbial
        if k_set.intersection(m_set): 

            # get compound from reactant side of equation
            for reactant in equation[0]: 

                # get compound from product side of equation
                for product in equation[1]: 

                    # ensure both product and reactant are nodes
                    if reactant in compounds and product in compounds: 
                        edges_dict['compound1'].append(reactant)
                        edges_dict['compound2'].append(product)
                        edges_dict['reaction'].append(rxn)

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
        else: 
            # a dict of compound associated with food KOs
            for compound in equation[0] + equation[1]:
                if compound not in comp_food:
                    comp_food[compound] = set()  # initialize empty set
                comp_food[compound].update(k_set)  # directly add new KOs

    edges_df = pd.DataFrame(edges_dict)
    edges_df.to_csv(o_edges, index=False)

    # create node information dictionaries 
    origin_dict = {}
    food_item_dict = {}
    food_freq_dict = {}

    # print some summary info
    print(f"There are {(mapper_df['origin'] == 'microbe').sum()} compounds assocaited with microbes")
    print(f"There are {(mapper_df['origin'] == 'food').sum()} compounds assocaited with food items")
    print(f"There are {(mapper_df['origin'] == 'both').sum()} compounds assocaited with both microbes and food items")

    # for each detected compound 
    for compound in compound_ids: 

        # get origin for compound
        origin = mapper_df.loc[mapper_df['id']==compound, 'origin'].values[0]
        origin_dict[compound] = origin

        # if compound is associated with a food 
        if origin_dict[compound] != 'microbe' and compound in comp_food:

            # get associated KOs 
            assoc_kos = comp_food[compound]

            comp_df = food_meta[food_meta['kegg_id'].isin(assoc_kos)]

            # get a list of unique food items 
            foods = set(comp_df['ScientificName'])
            food_item_dict[compound] = foods

            # sum frequency of consumption and make sure there is no value over 100 
            if n_weights: # if weights are being used 
                freq_sub = comp_df.drop_duplicates(subset=['ScientificName'])
                freq = freq_sub['food_frequency'].sum()
                if freq > 100:
                    raise ValueError(f"The frequency of consumption is too high: {freq}")
                food_freq_dict[compound] = int(freq)
            else:
                food_freq_dict[compound] = pd.NA # default no weights is NA 
        else: # microbe only associated compounds 
            food_item_dict[compound] = pd.NA
            food_freq_dict[compound] = pd.NA


    # Convert each dict to a DataFrame
    df_origin = pd.DataFrame(origin_dict.items(), columns=['c_id', 'origin'])
    df_food = pd.DataFrame(food_item_dict.items(), columns=['c_id', 'assoc_food'])
    df_freq = pd.DataFrame(food_freq_dict.items(), columns=['c_id', 'freq'])

    nodes_df = df_origin.merge(df_food, on='c_id').merge(df_freq, on='c_id')
    nodes_df.to_csv(o_nodes, index=False)
    return nodes_df, edges_df
