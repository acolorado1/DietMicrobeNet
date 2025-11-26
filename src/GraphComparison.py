import pandas as pd 
import ast 
import numpy as np
import matplotlib.pyplot as plt 
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
import argparse 
import os 

def get_graphs(paths:list, names:list): 
    """reads in all graph results CSVs and creates a dict with graphs associated with a name 

    Args:
        paths (list): list of file paths 
        names (list): list of names 

    Raises:
        ValueError: paths list and names list must be of same length
        ValueError: names must not be duplicated 

    Returns:
        dict: keys are names, values are pandas dataframes  
    """
    # make sure number of names and paths are equal 
    if len(names) != len(paths): 
        raise ValueError('Number of names not equal to the number of paths provided.')

    # Put graph results files is a list 
    graph_dict = {}
    for index in range(len(paths)):
        
        name = names[index]

        # ensure no duplicated names
        if name in graph_dict:
            raise ValueError('Name is duplicated, require all unique names')
        df = pd.read_csv(paths[index])
        graph_dict[name] = df

    return graph_dict

def subset_graphs(graph_dict:dict):
    """subset the overall graph dataframe into pattern specific results 

    Args:
        graph_dict (dict): graph_results CSV as pandas dataframe 

    Returns:
        dict, dict, dict: three dictionaries where keys are names and values are pandas dataframes 
    """
    #initialize dicts for each pattern type 
    food_microbe = {}
    food_both = {}
    both_both = {}

    # separate each graph into patterns 
    for name, df in graph_dict.items(): 
        food_microbe[name] = df[(df['compound1_origin']=='food') & (df['compound2_origin']=='microbe')]
        food_both[name] = df[(df['compound1_origin']=='food') & (df['compound2_origin']=='both')]
        both_both[name] = df[(df['compound1_origin']=='both') & (df['compound2_origin']=='both')]
    
    # return 3 dicts corresponding to each pattern detected 
    return food_microbe, food_both, both_both 

def get_kos(graph_dict:dict):
    """extract KOs from pandas dataframes representing the graphs 

    Args:
        graph_dict (dict): output of get_graphs()

    Returns:
        dict: keys are names, values list of KOs no duplicates 
    """
    # get kos for each graph
    kos_dict = {}
    for name, df in graph_dict.items(): 
        ko_column = df['KOs'].to_list() # format: ["['KOXXXX']"]

        # Parse string to list of lists 
        parsed = [ast.literal_eval(x) for x in ko_column] # format: [['KOXXXX']]

        # unnest the lists 
        unnested = [item for sublist in parsed for item in sublist] # format: ['KOXXXX']

        # remove duplicates 
        ko_set = list(set(unnested)) 

        #sort list 
        ko_set_sort = sorted(ko_set) 

        # append to dict 
        kos_dict[name] = ko_set_sort

    return kos_dict

def jaccard (a:list, b:list): 
    """calculate jaccard similarity score between two lists 

    Args:
        a (list): list of KOs
        b (list): list of KOs 

    Returns:
        float: intersection over union between two lists 
    """
    set_a = set(a)
    set_b = set(b)
    return len(set_a & set_b) / len(set_a | set_b)


def calculate_similarity_matrix(pattern_kos:dict): 
    """for each pattern dictionary create a matrix of jaccard similarity values 

    Args:
        pattern_kos (dict): keys are names, values are list of KOs 

    Returns:
        list: list of list representing a matrix 
    """
    labels = list(pattern_kos.keys())
    n = len(labels)
    matrix = np.zeros((n,n))

    for i, k1 in enumerate(labels):
        for j, k2 in enumerate(labels):
            matrix[i, j] = jaccard(pattern_kos[k1], pattern_kos[k2])

    return matrix, labels

def cluster_matrix(matrix:list, labels:list): 
    """from calculated matrix order using hierarchical clustering 

    Args:
        matrix (list): list of list representing a matrix 
        labels (list): list of names used as labels for plotting  

    Returns:
        list, list, list: list of lists representing an ordered matrix, 
                          list of ordered labels, 
                          linkage matrix returned by SciPy
    """
    distance_matrix = 1-matrix
    Z = linkage(distance_matrix, method='average') # hierarchical clustering 
    leaf_order = leaves_list(Z) # get leaf order 

    # reorder matrix and labels 
    clustered_matrix = matrix[leaf_order][:, leaf_order]
    clustered_labels = [labels[i] for i in leaf_order]

    return clustered_matrix, clustered_labels, Z 

def plotting(pattern_dict:dict, pattern_name:str, output:str):
    """given a dict containing names as keys and list of KOs as values 
    this script creates a similarity matrix and a dendrogram which are saved 

    Args:
        pattern_dict (dict): keys are names and values are list of KOs 
        pattern_name (str): name of the pattern (e.g., Food to Microbe)
        output (str): directory where output PNG files will go 
    """
    matrix, labels = calculate_similarity_matrix(pattern_dict)
    ordered_matrix, ordered_labels, Z = cluster_matrix(matrix=matrix, labels=labels)
    
    # Heatmap
    n = len(ordered_labels)
    plt.figure(figsize=(10, 8))
    plt.imshow(ordered_matrix, cmap='viridis')
    plt.xticks(range(n), ordered_labels, rotation=90)
    plt.yticks(range(n), ordered_labels)
    plt.colorbar(label="Jaccard Similarity")
    plt.title("Clustered Jaccard Similarity Heatmap: " + pattern_name)
    plt.tight_layout()
    plt.savefig(output + pattern_name.replace(" ", "") + 'GraphComparisons_Heatmap.png')
    plt.close()

    # Dendrogram
    plt.figure(figsize=(8, 4))
    dendrogram(Z, labels=ordered_labels, leaf_rotation=90)
    plt.title("Hierarchical Clustering Dendrogram: " + pattern_name)
    plt.tight_layout()
    plt.savefig(output + pattern_name.replace(" ", "") + 'GraphComparisons_Dendrogram.png')
    plt.close()

def summary(pattern_dict:dict, pattern_name:str, output:str):
    """create a summary text file of the graph comparisons 

    Args:
        pattern_dict (dict): values are names (e.g., sample IDs), and keys are list of KOs 
        pattern_name (str): identifies which pattern is being summarised 
        output (str): directory were output txt file will go 
    """
    # find intersection 
    intersection_set = set()
    for ko_list in pattern_dict.values():
        intersection_set.intersection(set(ko_list))

    unique_dict = {
        name: list(set(ko_list) - intersection_set)
        for name, ko_list in pattern_dict.items()
    }

    directory = output + pattern_name.replace(" ", "") + "GraphComparisons_Summary.txt"

    with open(directory, 'w') as file_object: 
        file_object.write(f'##### SUMMARY FOR PATTERN: {pattern_name} #####\n')
        file_object.write(f'Number of compounds shared: {len(list(intersection_set))}\n')

        for name, unique_list in unique_dict.items():
            file_object.write(f'Unique compounds to {name}: {len(unique_list)}\n')

        file_object.write(f'\n\n\n ###### LISTS OF KOs ######\n\n\n')
        file_object.write(f'Intersection KOs: {list(intersection_set)}\n\n')

        for name, unique_list in unique_dict.items(): 
            file_object.write(f'Unique KOs for {name}: {unique_list}\n\n')

    file_object.close()

# --- Argument parser ---
def main():
    parser = argparse.ArgumentParser(description="Compare graph results across samples using KOs and Jaccard similarity.")
    parser.add_argument("-p", "--paths", nargs='+', required=True, help="Paths to CSV graph results files")
    parser.add_argument("-n", "--names", nargs='+', required=True, help="Names corresponding to each CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output directory for plots and summary files")
    args = parser.parse_args()

    graphs_dict = get_graphs(paths=args.paths, names=args.names)
    food_microbe_dict, food_both_dict, both_both_dict = subset_graphs(graph_dict=graphs_dict)
    food_microbe_kos = get_kos(food_microbe_dict)
    food_both_kos = get_kos(food_both_dict)
    both_both_kos = get_kos(both_both_dict)

    patterns = [food_microbe_kos, food_both_kos, both_both_kos]
    pattern_names = ['Food to Microbe', 'Food to Both', 'Both to Both']

    for index in range(len(patterns)): 
        plotting(patterns[index], pattern_names[index], args.output)
        summary(patterns[index], pattern_names[index], args.output)

if __name__ == "__main__":
    main()
