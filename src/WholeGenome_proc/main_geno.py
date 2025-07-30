import argparse
import geno_nodes_edges as gne

def parse_args():
    parser = argparse.ArgumentParser(description="Create node and edge CSVs from genome predictions and AMON outputs.")

    parser.add_argument('--f_meta', type=str, required=True,
                        help="File path to CSV created from comp_KEGG.py (food info)")
    parser.add_argument('--m_meta', type=str, required=True,
                        help="File path to CSV containing KO, taxonomic, and abundance info")
    parser.add_argument('--mapper', type=str, required=True,
                        help="Mapper file path generated when running AMON")
    parser.add_argument('--co_json', type=str, required=True,
                        help="JSON file path generated on compound information from AMON")
    parser.add_argument('--rn_json', type=str, required=True,
                        help="JSON file path generated on reaction information from AMON")
    parser.add_argument('--n_weights', action='store_true',
                        help="If provided, nodes will have frequencies associated with them")
    parser.add_argument('--e_weights', action='store_true',
                        help="If provided, edges will have frequencies associated with them")
    parser.add_argument('--o_nodes', type=str, required=True,
                        help="Output file path for nodes CSV")
    parser.add_argument('--o_edges', type=str, required=True,
                        help="Output file path for edges CSV")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    gne.node_edge_dfs(
        f_meta=args.f_meta,
        m_meta=args.m_meta,
        mapper=args.mapper,
        co_json=args.co_json,
        rn_json=args.rn_json,
        n_weights=args.n_weights,
        e_weights=args.e_weights,
        o_nodes=args.o_nodes,
        o_edges=args.o_edges
    )

