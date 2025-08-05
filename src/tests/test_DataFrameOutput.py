import unittest
from pathlib import Path
from Metabolome_proc import metab_nodes_edges as Mne
from WholeGenome_proc import geno_nodes_edges as Gne
import os

class MyTestCase(unittest.TestCase): 

    # make sure column names are as expected 
    def test_MetabDataframes(self):
        project_root = Path(__file__).resolve().parents[2] # this takes us to DietMicrobeNet 
        microbe_comps_path = project_root / "Data" / "metab_test" / "AMON_output" / "co_dict.json"
        food_comps_path = project_root / "Data" / "test" / "metab_outputs" / "compound_info" / "foodb_meta_freq.csv"
        rn_path = project_root / "Data" / "metab_test" / "AMON_output" / "rn_dict.json"
        ko_meta_path = project_root / "Data" / "sample" / "ko_taxonomy_abundance.csv"

        nodes = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=True)

        edges = Mne.edge_df(rn_json=rn_path, 
                                      ko_meta=ko_meta_path, 
                                      nodes_df=nodes, 
                                      e_weights=True)
        
        # list of expected columns 
        node_columns = ['c_id', 'origin', 'assoc_food', 'freq']
        edge_columns = ['compound1','compound2','reaction','KOs','organisms','abundance_RPKs']

        # check column values 
        self.assertEqual(list(nodes.columns), node_columns)
        self.assertEqual(list(edges.columns), edge_columns)


    # test genome node and edge weights 
    def test_GenoNodeAndEdgeWeights(self): 
        project_root = Path(__file__).resolve().parents[2] # this takes us to DietMicrobeNet 
        f_meta_path = project_root / "Data" / "test" / "org_KO_adj" / "orgs_KOs.csv"
        m_meta_path = project_root / "Data" / "sample" / "ko_taxonomy_abundance.csv"
        mapper_path = project_root / "Data" / "geno_test" / "AMON_output" / "kegg_mapper.tsv"
        co_path = project_root / "Data" / "geno_test" / "AMON_output" / "co_dict.json"
        rn_path = project_root / "Data" / "geno_test" / "AMON_output" / "rn_dict.json"

        nodes, edges = Gne.node_edge_dfs(f_meta=f_meta_path, 
                                             m_meta=m_meta_path,
                                             mapper=mapper_path, 
                                             co_json=co_path,
                                             rn_json=rn_path,
                                             n_weights=True,
                                             e_weights=True,
                                             o_nodes='nodes.csv',
                                             o_edges='edges.csv')
        # list of expected columns 
        node_columns = ['c_id', 'origin', 'assoc_food', 'freq']
        edge_columns = ['compound1','compound2','reaction','KOs','organisms','abundance_RPKs']

        # check column values 
        self.assertEqual(list(nodes.columns), node_columns)
        self.assertEqual(list(edges.columns), edge_columns)

        # remove the files created by this test 
        os.remove(project_root / "nodes.csv")
        os.remove(project_root / "edges.csv")

if __name__ == "__main__":
    unittest.main()