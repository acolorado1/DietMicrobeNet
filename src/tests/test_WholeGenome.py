import unittest
import pandas as pd
from WholeGenome_proc import nodes_edges as ne 

# create dummy data 
food_meta_df = pd.DataFrame({
    'kegg_id': ['KOOOO1', 'K00002', 'K00002', 'K00003'], 
    'ScientificName': ['Apple', 'Apple', 'Cow', 'Cow'],
    'food_frequency': [60, 60, 40, 40]
})

microbe_meta = pd.DataFrame({
    'KO': ['K00001', 'K00001', 'K00004'],
    'taxonomy': ['org1', 'org2', 'org3'],
    'Abundance_RPKs': [5, 10, 50]
})

mapper = pd.DataFrame({
    'id': ['K00001', 'K00002', 'K00003', 'K00004', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9'],
    'amon_origin': ['green', 'yellow', 'yellow', 'blue', 'green', 'green', 'green', 'yellow', 'green', 'green', 'green', 'blue', 'blue'],
    'origin':      ['both', 'food', 'food', 'microbe', 'both', 'both', 'both', 'food', 'both', 'both', 'both', 'microbe', 'microbe']
})

class MyTestCase(unittest.TestCase): 
    def test_org_abundance(self): 
        abundance, orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          abundance_column='Abundance_RPKs')
        expected_abundance = {'K00001': 15, 'K00004': 50}
        expected_orgs = {'K00001': 'org1, org2', 'K00004': 'org3'}  
        
        self.assertEqual(abundance, expected_abundance)
        self.assertEqual(orgs, expected_orgs)
    
    def test_info_dicts(self):
        comp_kos, rn_kos, rn_eq = ne.get_info_dicts(rxns=rn_dict)
        expected_compkos_c1 = ['K00001', 'K00003']
        expected_rnkos_rn4 = ['K00003', 'K00004']
        expected_rn2 = [['C2', 'C3'], ['C4', 'C5']]

        self.assertCountEqual(comp_kos['C1'], expected_compkos_c1)
        self.assertCountEqual(rn_kos['rn4'], expected_rnkos_rn4)
        self.assertCountEqual(rn_eq['rn2'], expected_rn2)

    def test_get_organisms_and_abundance(self): 
        kos = ['K00003', 'K00004']
        abundance, orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          abundance_column='Abundance_RPKs')

        o = ne.get_organisms(kos, orgs)
        a = ne.get_abundance(kos, abundance)
        
        expected_o = ['org3']
        expected_a = 50 

        self.assertCountEqual(o , expected_o)
        self.assertEqual(a, expected_a)
    
    def test_edge_creation(self): 
        abundance, orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          abundance_column='Abundance_RPKs')
        comp_kos, rn_kos, rn_eq = ne.get_info_dicts(rxns=rn_dict)

        edge_df = ne.build_edges_df(mapper=mapper, 
                                     rxn_kos=rn_kos, 
                                     rxn_equation=rn_eq,  
                                     orgs=True, 
                                     e_weights=True, 
                                     ko_organisms=orgs, 
                                     ko_abundance=abundance)

        # test how many edges would be made, expect 5
        self.assertEqual(len(edge_df), 5) 
        
        # test dataframe output 
        edge_columns = ['compound1','compound2','reaction','KOs','organisms','abundance']
        self.assertEqual(list(edge_df.columns), edge_columns)

        # test that edge weights become NAs 
        edge_wo_weights = ne.build_edges_df(mapper=mapper, 
                                     rxn_kos=rn_kos, 
                                     rxn_equation=rn_eq, 
                                     orgs=True, 
                                     e_weights=False, 
                                     ko_organisms=orgs, 
                                     ko_abundance=abundance)
        self.assertTrue(edge_wo_weights['abundance'].isnull().all())

        # test that organisms become NAs 
        edge_wo_orgs = ne.build_edges_df(mapper=mapper, 
                                     rxn_kos=rn_kos, 
                                     rxn_equation=rn_eq, 
                                     orgs=False, 
                                     e_weights=True, 
                                     ko_organisms=orgs, 
                                     ko_abundance=abundance)
        self.assertTrue(edge_wo_orgs['organisms'].isnull().all())

    def test_node_creation(self):
        abundance, orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          abundance_column='Abundance_RPKs')
        comp_kos, rn_kos, rn_eq = ne.get_info_dicts(rxns=rn_dict) 

        node_df = ne.build_nodes_df(mapper=mapper, 
                                     food_meta=food_meta_df,
                                     compound_kos=comp_kos,
                                     frequency=True)

        # test how many nodes would be made, expect 9
        self.assertEqual(len(node_df), 9) 

        # test that all frequencies are either 60, 40, 100 or NA
        self.assertTrue(node_df['freq'].isin([40, 60, 100, pd.NA]).all())

        # make sure that all microbial compounds have no associated food 
        associations = node_df[node_df['origin']=='microbe']
        check = associations['food_assoc'].isna().all() 
        self.assertTrue(check)

        # make sure all microbial compound have no associated freq 
        associations = node_df[node_df['origin']=='microbe']
        check = associations['freq'].isna().all() 
        self.assertTrue(check)

        # make sure all food associated compounds have food assoc 
        associations = node_df[node_df['origin']!='microbe']
        check = associations['food_assoc'].isna().all() 
        self.assertFalse(check)

        # make sure all microbial compound have no associated freq 
        associations = node_df[node_df['origin']!='microbe']
        check = associations['freq'].isna().all() 
        self.assertFalse(check)

        # make sure weights are na when frequency = False
        node_df = ne.build_nodes_df(mapper=mapper, 
                                     food_meta=food_meta_df,
                                     compound_kos=comp_kos,
                                     frequency=False)
        check = node_df['freq'].isna().all() 
        self.assertTrue(check)

if __name__ == "__main__":
    unittest.main()
