import unittest
from pathlib import Path
from Metabolome_proc import metab_nodes_edges as Mne

project_root = Path(__file__).resolve().parents[2] # this takes us to DietMicrobeNet 
microbe_comps_path = project_root / "Data" / "metab_test" / "AMON_output" / "co_dict.json"
food_comps_path = project_root / "Data" / "test" / "metab_outputs" / "compound_info" / "foodb_meta_freq.csv"
rn_path = project_root / "Data" / "metab_test" / "AMON_output" / "rn_dict.json"
ko_meta_path = project_root / "Data" / "sample" / "ko_taxonomy_abundance.csv"

class MyTestCase(unittest.TestCase): 

    # make sure column names are as expected 
    def test_MetabDataframes(self):

        nodes = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=True)

        edges = Mne.edge_df(rn_json=rn_path, 
                                      ko_meta=ko_meta_path, 
                                      nodes_df=nodes, 
                                      e_weights=True, 
                                      abundance='Abundance_RPKs',
                                      orgs=True)
        
        # list of expected columns 
        node_columns = ['c_id', 'origin', 'assoc_food', 'freq']
        edge_columns = ['compound1','compound2','reaction','KOs','organisms','abundance']

        # check column values 
        self.assertEqual(list(nodes.columns), node_columns)
        self.assertEqual(list(edges.columns), edge_columns)

    # test with and without node weights  
    def test_MetabNodeWeights(self): 

        nodes_w_weights = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=True)

        nodes_wo_wights = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=False)
    
        # check column values 
        self.assertFalse(nodes_w_weights['freq'].isna().all()) 
        self.assertTrue(nodes_wo_wights['freq'].isna().all())

    # test with and without edge weights 
    def test_MetabEdgeWeights(self):

        nodes = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=True)

        edges_w_weights = Mne.edge_df(rn_json=rn_path, 
                                      ko_meta=ko_meta_path, 
                                      nodes_df=nodes, 
                                      e_weights=True, 
                                      abundance='Abundance_RPKs',
                                      orgs=True)
        edges_wo_weights = Mne.edge_df(rn_json=rn_path, 
                                      ko_meta=ko_meta_path, 
                                      nodes_df=nodes, 
                                      e_weights=False, 
                                      abundance='Abundance_RPKs',
                                      orgs=True)
        
        # check column values 
        self.assertFalse(edges_w_weights['abundance'].isna().all()) 
        self.assertTrue(edges_wo_weights['abundance'].isna().all())

    # test metabolome node origins
    def test_MetabNodeOrigins(self):   

        nodes = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=True)
        
        # list of only possible origins 
        possible_origins = ['both', 'microbe', 'food']

        # get origins column
        origins = nodes['origin']

        self.assertTrue(origins.isin(possible_origins).all())

if __name__ == "__main__":
    unittest.main()