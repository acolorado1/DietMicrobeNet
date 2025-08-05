import unittest
from pathlib import Path
from Metabolome_proc import metab_nodes_edges as Mne
from WholeGenome_proc import geno_nodes_edges as Gne
import os

class MyTestCase(unittest.TestCase): 

    # test metabolome node origins
    def test_MetabNodeOrigins(self):   # sourcery skip: class-extract-method
        project_root = Path(__file__).resolve().parents[2] # this takes us to DietMicrobeNet 
        microbe_comps_path = project_root / "Data" / "metab_test" / "AMON_output" / "co_dict.json"
        food_comps_path = project_root / "Data" / "test" / "metab_outputs" / "compound_info" / "foodb_meta_freq.csv"

        nodes = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=True)
        
        # list of only possible origins 
        possible_origins = ['both', 'microbe', 'food']

        # get origins column
        origins = nodes['origin']

        self.assertTrue(origins.isin(possible_origins).all())



    # test genome node origins 
    def test_GenoNodeOrigins(self): 
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
        
        # list of only possible origins 
        possible_origins = ['both', 'microbe', 'food']

        # get origins column
        origins = nodes['origin']

        self.assertTrue(origins.isin(possible_origins).all())

        # remove the files created by this test 
        os.remove(project_root / "nodes.csv")
        os.remove(project_root / "edges.csv")

if __name__ == "__main__":
    unittest.main()