# this code will run the edges and nodes dataframe creation from terminal 

import metab_nodes_edges as mne 
import argparse as arg

def main():
    parser = arg.ArgumentParser(
        description='Output two dataframes as CSVs of node and edge information'
    )
    parser.add_argument(
        '-m',
        '--microbe_info',
        type=str, 
        required=True,
        help='file path to AMON output of microbial compounds, KOs, and reactions'
    )
    parser.add_argument(
        '-f', 
        '--food_info', 
        type=str,
        required=True, 
        help='file path to food compound information taken from comp_FoodDB.R'
    )
    parser.add_argument(
        '-n_w', 
        '--node_weights', 
        action='store_true', 
        required=True, 
        help='True or False wether node weights (food frequencies) will be used'
    )
    parser.add_argument( 
        '-r', 
        '--reactions', 
        type=str, 
        required=True, 
        help='file path to rn_dict.json AMON output'
    )
    parser.add_argument(
        '-meta', 
        '--microbiome_metadata', 
        type=str, 
        required=True, 
        help='filepath to CSV containing KOs, taxonomy, and abundance information'
    )
    parser.add_argument(
        '-e_w',
        '--edge_weights', 
        action='store_true', 
        required=True, 
        help='True or False wether edge weights (abundance measures) will be used'
    )
    parser.add_argument(
        '-n_o', 
        '--nodes_output', 
        type=str, 
        required=True, 
        help='directory where CSV of node information should go'
    )
    parser.add_argument(
        '-e_o', 
        '--edges_output', 
        type=str, 
        required=True, 
        help='directory where CSV of edge information should go'
    )

    args = parser.parse_args()

    # call nodes and edges functions 
    node_df = mne.node_df(microbome_compounds=args.microbe_info, 
                           food_compounds=args.food_info, 
                           n_weights=args.node_weights)
    edge_df = mne.edge_df(rn_json=args.reactions, 
                           ko_meta=args.microbiome_metadata, 
                           nodes_df=node_df, 
                           e_weights=args.edge_weights)
    # create CSVs 
    node_df.to_csv(args.nodes_output, index=False)
    edge_df.to_csv(args.edges_output, index=False)

    exit()

if __name__ == "__main__":
    main()