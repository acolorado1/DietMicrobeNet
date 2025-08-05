import unittest
from pathlib import Path
from Metabolome_proc import metab_nodes_edges as Mne
from WholeGenome_proc import geno_nodes_edges as Gne
import os

class MyTestCase(unittest.TestCase): 

    # test metabolome node and edge weights 
    def test_MetabNodeWeights(self): 
        project_root = Path(__file__).resolve().parents[2] # this takes us to DietMicrobeNet 
        microbe_comps_path = project_root / "Data" / "metab_test" / "AMON_output" / "co_dict.json"
        food_comps_path = project_root / "Data" / "test" / "metab_outputs" / "compound_info" / "foodb_meta_freq.csv"

        nodes_w_weights = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=True)

        nodes_wo_wights = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=False)
    
        # check column values 
        self.assertFalse(nodes_w_weights['freq'].isna().all()) 
        self.assertTrue(nodes_wo_wights['freq'].isna().all())

    def test_MetabEdgeWeights(self):
        project_root = Path(__file__).resolve().parents[2] # this takes us to DietMicrobeNet 
        microbe_comps_path = project_root / "Data" / "metab_test" / "AMON_output" / "co_dict.json"
        food_comps_path = project_root / "Data" / "test" / "metab_outputs" / "compound_info" / "foodb_meta_freq.csv"
        rn_path = project_root / "Data" / "metab_test" / "AMON_output" / "rn_dict.json"
        ko_meta_path = project_root / "Data" / "sample" / "ko_taxonomy_abundance.csv"

        nodes = Mne.node_df(microbome_compounds=microbe_comps_path, 
                            food_compounds=food_comps_path,
                            n_weights=True)

        edges_w_weights = Mne.edge_df(rn_json=rn_path, 
                                      ko_meta=ko_meta_path, 
                                      nodes_df=nodes, 
                                      e_weights=True)
        edges_wo_weights = Mne.edge_df(rn_json=rn_path, 
                                      ko_meta=ko_meta_path, 
                                      nodes_df=nodes, 
                                      e_weights=False)
        
        # check column values 
        self.assertFalse(edges_w_weights['abundance_RPKs'].isna().all()) 
        self.assertTrue(edges_wo_weights['abundance_RPKs'].isna().all())


    # test genome node and edge weights 
    def test_GenoNodeAndEdgeWeights(self): 
        project_root = Path(__file__).resolve().parents[2] # this takes us to DietMicrobeNet 
        f_meta_path = project_root / "Data" / "test" / "org_KO_adj" / "orgs_KOs.csv"
        m_meta_path = project_root / "Data" / "sample" / "ko_taxonomy_abundance.csv"
        mapper_path = project_root / "Data" / "geno_test" / "AMON_output" / "kegg_mapper.tsv"
        co_path = project_root / "Data" / "geno_test" / "AMON_output" / "co_dict.json"
        rn_path = project_root / "Data" / "geno_test" / "AMON_output" / "rn_dict.json"

        w_nodes, w_edges = Gne.node_edge_dfs(f_meta=f_meta_path, 
                                             m_meta=m_meta_path,
                                             mapper=mapper_path, 
                                             co_json=co_path,
                                             rn_json=rn_path,
                                             n_weights=True,
                                             e_weights=True,
                                             o_nodes='nodes.csv',
                                             o_edges='edges.csv')
        
        wo_nodes, wo_edges = Gne.node_edge_dfs(f_meta=f_meta_path, 
                                             m_meta=m_meta_path,
                                             mapper=mapper_path, 
                                             co_json=co_path,
                                             rn_json=rn_path,
                                             n_weights=False,
                                             e_weights=False,
                                             o_nodes='nodes.csv',
                                             o_edges='edges.csv')
        
        # check nodes 
        self.assertFalse(w_nodes['freq'].isna().all()) 
        self.assertTrue(wo_nodes['freq'].isna().all())
        
        # check edges 
        self.assertFalse(w_edges['abundance_RPKs'].isna().all()) 
        self.assertTrue(wo_edges['abundance_RPKs'].isna().all())

        # remove the files created by this test 
        os.remove(project_root / "nodes.csv")
        os.remove(project_root / "edges.csv")

if __name__ == "__main__":
    unittest.main()