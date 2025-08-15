import argparse
import nodes_edges as ne
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import sys
import warnings

def parse_args():
    parser = argparse.ArgumentParser(description="Create node and edge CSVs from genome predictions and AMON outputs.")

    parser.add_argument('--f_meta', type=str, required=True,
                        help="File path to CSV created from comp_KEGG.py (food info)")
    parser.add_argument('--m_meta', type=str, required=True,
                        help="File path to CSV containing KO, taxonomic, and abundance info")
    parser.add_argument('--mapper', type=str, required=True,
                        help="Mapper file path generated when running AMON")
    parser.add_argument('--rn_json', type=str, required=True,
                        help="JSON file path generated on reaction information from AMON")
    parser.add_argument('--n_weights', action='store_true',
                        help="If provided, nodes will have frequencies associated with them")
    parser.add_argument('--e_weights', action='store_true',
                        help="If provided, edges will have frequencies associated with them")
    parser.add_argument('--a', type=str, required=True, default=None,
                        help='name of the column where abundance information is located')
    parser.add_argument('--org', action='store_true',
                        help='If provided, will add organism information to edge dataframe')
    parser.add_argument('--o', type=str, required=True,
                        help="Output file path for nodes CSV")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # check if output directory exists and create one if it does not 
    if not os.path.exists(args.o):
        os.makedirs(args.o)
        print(f'Directory {args.o} created successfully')
    else:
        print('Directory already exists')
        sys.exit()

    # read in all files 
    food_meta, m_meta_clean, mapper_df, rxns = ne.data_read_in(f_meta=args.f_meta, 
                                                                m_meta=args.m_meta, 
                                                                mapper=args.mapper, 
                                                                rn_json=args.rn_json, 
                                                                abundance_column=args.a, 
                                                                e_weights=args.e_weights, 
                                                                orgs=args.org)
    
    # get abundance and organism dicts
    ko_abundance, ko_orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=m_meta_clean,
                                                              abundance_column=args.a, 
                                                              e_weights=args.e_weights, 
                                                              orgs=args.org)
    
    # get all needed info from reaction json 
    compound_kos, rxn_kos, rxn_equation = ne.get_info_dicts(rxns=rxns)

    # create edges pandas dataframe 
    edges_df = ne.build_edges_df(mapper=mapper_df,
                                  rxn_kos=rxn_kos, 
                                  rxn_equation=rxn_equation,
                                  orgs=args.org, 
                                  e_weights=args.e_weights,
                                  ko_organisms=ko_orgs, 
                                  ko_abundance=ko_abundance)
    
    if len(edges_df) == 0:
        warnings.warn('There are no edges found, please check files to make sure they are correct ones')
    
    # create nodes dataframe 
    nodes_df = ne.build_nodes_df(mapper=mapper_df,
                                  food_meta=food_meta, 
                                  compound_kos=compound_kos, 
                                  frequency=args.n_weights)
    if len(nodes_df) == 0:
        warnings.warn('There are no nodes found, please check files to make sure they are correct ones')
    
    # have summary print
    ne.summarize_res(compound_kos=compound_kos, 
                      rxn_kos=rxn_kos,
                      nodes_df=nodes_df, 
                      edges_df=edges_df) 
    
    # create distribution plots of food frequency and abundance 
    if args.n_weights:
        plt.figure(figsize=(8,6))
        sns.histplot(nodes_df['freq'], kde=True)
        plt.title('Distribution of Food Frequencies Associated with Nodes')
        plt.xlabel('Food Frequency Value')
        plt.ylabel('Frequency/Density')
        plt.savefig(args.o + 'WG_FoodFrequencyDistribution.png')
        plt.close()

    if args.e_weights:
        plt.figure(figsize=(8,6))
        sns.histplot(edges_df['abundance'], kde=True)
        plt.title('Distribution of KO Abundance Associated with Edges')
        plt.xlabel('Abundance')
        plt.ylabel('Frequency/Density')
        plt.savefig(args.o + 'WG_AbundanceDistribution.png')
        plt.close()
    
    # output as csvs 
    edges_df.to_csv(args.o + 'WG_edges_df.csv', index=False)
    nodes_df.to_csv(args.o + 'WG_nodes_df.csv', index=False)

