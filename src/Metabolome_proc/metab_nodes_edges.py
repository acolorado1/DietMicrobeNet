import json
import pandas as pd

def node_df(microbome_compounds:str, food_compounds:str):
    """create node dataframe used as input for graph visualization 

    Args:
        microbome_compounds (str): path to compound json file from AMON output 
        food_compounds (str): path to CSV output from comp_FoodDB.py 
    
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

            # sum frequency of consumption and make sure there is no value over 100 
            freq_sub = comp_df.drop_duplicates(subset=['name'])
            freq = freq_sub['food_frequency'].sum()
            if freq > 100:
                raise ValueError(f"The frequency of consumption is too high: {freq}")

            
            # add to dictionaries 
            food_item_dict[compound] = foods
            food_freq_dict[compound] = int(freq)

    # Convert each dict to a DataFrame
    df_origin = pd.DataFrame(origin_dict.items(), columns=['c_id', 'origin'])
    df_food = pd.DataFrame(food_item_dict.items(), columns=['c_id', 'assoc_food'])
    df_freq = pd.DataFrame(food_freq_dict.items(), columns=['c_id', 'freq'])

    # Merge them on 'c_id'
    df_combined = df_origin.merge(df_food, on='c_id').merge(df_freq, on='c_id')


    print(df_combined.head())
 

# Files needed for node creation 
#m = '/Users/burkhang/Code_Projs/DietMicrobeNet/Data/test/AMON_metab_output/co_dict.json' #list of microbe compounds 
#f = '/Users/burkhang/Code_Projs/DietMicrobeNet/Data/test/compound_meta/foodb_meta_freq.csv' # list of food compounds 
#node_df(m, f)