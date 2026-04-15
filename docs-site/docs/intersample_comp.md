# Inter-Sample Comparison 

## Comparing N Graphs/Pattern Outputs 

After running all previous steps you will end up with a file `/graph/graph_results.csv` which contains the results of three different queries to find instances of microbial metabolism of dietary compounds. 

In order to find similarities and differences between the graphs/patterns we looked at the genes invovled in this metabolism and compared them using [Jaccard Similarity](https://mayurdhvajsinhjadeja.medium.com/jaccard-similarity-34e2c15fb524). This is vizualized in two ways, a heatmap of similarity scores, and a dendrogram to identify clusters using [SciPy's higherarchical clustering](https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html) algorithm. Additionally, a summary text file is written to show common genes between all graph for each pattern type and unique genes to each graph for each pattern type. The statistical test included in the summary is a [PERMANOVA](https://scikit.bio/docs/dev/generated/skbio.stats.distance.permanova.html) which performs 5,000 permutations and a seed of 5. 

## Running Comparison 

To get a list of optional and required arguments run **python GraphComparison.py -h**: 

```
python src/GraphComparison.py -h
usage: GraphComparison.py [-h] -m METADATA -p PATHS -n NAMES [-s] [-g GROUPS] -o OUTPUT [--ko_column KO_COLUMN]

Compare graph results across samples using KOs and Jaccard similarity.

options:
  -h, --help            show this help message and exit
  -m, --metadata METADATA
                        Metadata CSV containing file paths and names
  -p, --paths PATHS     Name of column containing file paths
  -n, --names NAMES     Name of column containing names of graphs (e.g., sampleID)
  -s, --stat_test       If statistical test for group comparison wanted include this parameter
  -g, --groups GROUPS   Names of columns for use in PERMANOVA, if multiple separate by a comma e.g., cohort,diet,location
  -o, --output OUTPUT   Output directory for plots and summary files
  --ko_column KO_COLUMN
                        Name of KOs column in graph CSVs (default: 'KOs')
```

NOTE 1: all `graph_results.csv` should be formatted the same unless there has been a change done by the user. Thus all the ko_columns column should use the default column name unless they have been changed.

NOTE 2: example metadata can be located in the `Data` file called **Example_GraphComparison_Metadata.csv**

Example usage w/o stats: 

```
python src/GraphComparison.py \ 
-m "path/to/metadata.csv" \
-p "paths_column_name" \
-n "names_column_name" \
-o "path/to/output/dir/" \
```

Example usage w/ stats: 

```
python src/GraphComparison.py \ 
-m "path/to/metadata.csv" \
-p "paths_column_name" \
-n "names_column_name" \
-s \
-g "groups_column_name" \
--o "path/to/output/dir/" \
```