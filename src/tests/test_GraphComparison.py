import unittest
import pandas as pd
import numpy as np
import GraphComparison as gc 

class MyTestCase(unittest.TestCase): 

    def setUp(self):
        sample1 = {'compound1_origin': ['food', 'food', 'both'],
                   'compound2_origin': ['microbe', 'both', 'both'], 
                   'KOs': ["['KO1']", "['KO2']", "['KO3']"]}
        sample2 = {'compound1_origin': ['food', 'food', 'both'], 
                   'compound2_origin': ['microbe', 'both', 'both'], 
                   'KOs': ["['KO1']", "['KO4']", "['KO5']"]}
        
        self.sample_dict = {
            'sample1': pd.DataFrame(sample1),
            'sample2': pd.DataFrame(sample2)
        }

    # Is it subsetting properly by pattern? 
    def test_subset(self): 
        fm, fb, bb = gc.subset_graphs(self.sample_dict)

        # Check rows in food_microbe
        self.assertEqual(fm['sample1'].shape[0], 1)
        self.assertEqual(fm['sample2'].shape[0], 1)
        
        # Check rows in food_both
        self.assertEqual(fb['sample1'].shape[0], 1)
        self.assertEqual(fb['sample2'].shape[0], 1)
        
        # Check rows in both_both
        self.assertEqual(bb['sample1'].shape[0], 1)
        self.assertEqual(bb['sample2'].shape[0], 1)

    # Does it grab the correct KOs and format?
    def test_get_kos(self): 
        fm, fb, bb = gc.subset_graphs(self.sample_dict)

        kos = gc.get_kos(fm)

        s1_value = kos['sample1']
        s2_value = kos['sample2']

        self.assertEqual(s1_value, ['KO1'])
        self.assertEqual(s2_value, ['KO1'])

    def test_jaccard(self): 
        a = ['KO1', 'KO2', 'KO3']
        b = ['KO1', 'KO4', 'KO5']

        score = gc.jaccard(a, b)
        expected = 1 / 5
        self.assertAlmostEqual(score, expected)

    def test_matrix_calc(self): 
        fm, fb, bb = gc.subset_graphs(self.sample_dict)

        ## check fm similarity matrix 
        kos = gc.get_kos(fm)
        matrix, labels = gc.calculate_similarity_matrix(kos)

        expected_matrix = np.array([
            [1.0, 1.0],
            [1.0, 1.0]
        ])

        np.testing.assert_array_equal(matrix, expected_matrix)

        # Also check that labels match input ordering
        self.assertEqual(labels, ['sample1', 'sample2'])

        ## check fb similarity matrix 
        kos = gc.get_kos(fb)
        matrix, labels = gc.calculate_similarity_matrix(kos)

        expected_matrix = np.array([
            [1.0, 0.0],
            [0.0, 1.0]
        ])

        np.testing.assert_array_equal(matrix, expected_matrix)

        # Also check that labels match input ordering
        self.assertEqual(labels, ['sample1', 'sample2'])

        ## check bb similarity matrix 
        kos = gc.get_kos(bb)
        matrix, labels = gc.calculate_similarity_matrix(kos)

        expected_matrix = np.array([
            [1.0, 0.0],
            [0.0, 1.0]
        ])

        np.testing.assert_array_equal(matrix, expected_matrix)

        # Also check that labels match input ordering
        self.assertEqual(labels, ['sample1', 'sample2'])


if __name__ == "__main__":
    unittest.main(argv=[''], exit=False)
