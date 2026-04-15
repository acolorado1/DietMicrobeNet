# Metabolome Comparison

## Comparing Pattern Results to a Provided Metabolome

To view which compounds identified in the patterns were also found in a metabolomics experiment you performed, two inputs are needed: 

### Inputs 

1. A `graph_results.csv` or the output of **python src/run_graph.py** in Step 4. 

2. CSV containing a list of KEGG compounds that were identified in the metabolome. An example of this file can be found in `Data/test_sample/metabolome.csv`.

### Running the script 

To get a list of optional and required arguments run **python src/RenderMetabolomeComparison.py -h**:

```
options:
  -h, --help            show this help message and exit
  --patterns PATTERNS   Path to the graph_results.csv
  --metabolome METABOLOME
                        Path to CSV file containing one column of KEGG compounds.
  --output OUTPUT       Path to HTML report file
```