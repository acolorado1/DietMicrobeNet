from __future__ import annotations
import argparse
import ast
import logging
import os
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from skbio import DistanceMatrix
from skbio.stats.distance import permanova

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def csv_to_inputs(metadata: str, paths_col: str, names_col: str, groups_col: str) -> Tuple[List[str], List[str], Dict[str, str]]:
    """Read metadata CSV and return paths, names, and optional group dict."""
    df = pd.read_csv(metadata)
    for col in (paths_col, names_col):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in metadata CSV ({metadata}). Available: {list(df.columns)}")

    paths_list = df[paths_col].astype(str).tolist()
    names_list = df[names_col].astype(str).tolist()

    if groups_col:
        if groups_col not in df.columns:
            raise ValueError(f"Groups column '{groups_col}' not found in metadata CSV ({metadata}).")
        # create dict name -> group
        group_dict = pd.Series(df[groups_col].astype(str).values, index=df[names_col].astype(str)).to_dict()
    else:
        group_dict = {}

    if len(paths_list) != len(names_list):
        raise ValueError("Number of paths and names differ in metadata CSV.")

    return paths_list, names_list, group_dict


def get_graphs(paths: List[str], names: List[str]) -> Dict[str, pd.DataFrame]:
    """Read each CSV path into a dataframe keyed by the corresponding name."""
    if len(paths) != len(names):
        raise ValueError("Number of names not equal to the number of paths provided.")

    graph_dict: Dict[str, pd.DataFrame] = {}
    for idx, p in enumerate(paths):
        name = names[idx]
        if name in graph_dict:
            raise ValueError(f"Duplicate name detected: '{name}'. Names must be unique.")
        pth = Path(p)
        if not pth.exists():
            raise FileNotFoundError(f"Graph CSV not found for sample '{name}': {p}")
        df = pd.read_csv(pth)
        graph_dict[name] = df

    return graph_dict


def subset_graphs(graph_dict: Dict[str, pd.DataFrame]) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """Split graph results into the three patterns used downstream.

    Validates required columns exist; raises informative errors if not.
    """
    required_cols = {"compound1_origin", "compound2_origin"}
    for name, df in graph_dict.items():
        if not required_cols.issubset(df.columns):
            raise ValueError(f"Missing required columns in graph dataframe for sample '{name}'. Required: {required_cols}. Found: {set(df.columns)}")

    food_microbe = {}
    food_both = {}
    both_both = {}

    for name, df in graph_dict.items():
        food_microbe[name] = df[(df["compound1_origin"] == "food") & (df["compound2_origin"] == "microbe")]
        food_both[name] = df[(df["compound1_origin"] == "food") & (df["compound2_origin"] == "both")]
        both_both[name] = df[(df["compound1_origin"] == "both") & (df["compound2_origin"] == "both")]

    return food_microbe, food_both, both_both


def _safe_literal_eval(item):
    """Try to parse a string representation of a list; return empty list on failure."""
    if item is None or (isinstance(item, float) and np.isnan(item)):
        return []
    # If it's already a list/tuple, return as list
    if isinstance(item, (list, tuple, set)):
        return list(item)
    # Strings: try literal_eval, fallback to simple parsing
    s = str(item).strip()
    if s == "":
        return []
    try:
        parsed = ast.literal_eval(s)
        if isinstance(parsed, (list, tuple, set)):
            return list(parsed)
        # if a single value (str) returned, wrap in list
        return [parsed]
    except Exception:
        # fallback: remove outer brackets/quotes and split on commas
        trimmed = s.strip("[](){} ")
        if trimmed == "":
            return []
        return [x.strip().strip("'\"") for x in trimmed.split(",") if x.strip()]


def get_kos(graph_dict: Dict[str, pd.DataFrame], ko_column_name: str = "KOs") -> Dict[str, List[str]]:
    """Extract KOs for each sample from the graph dataframes.

    Ensures unique, sorted KOs for determinism.
    """
    kos_dict: Dict[str, List[str]] = {}
    for name, df in graph_dict.items():
        if ko_column_name not in df.columns:
            logging.warning(f"'{ko_column_name}' column not found in dataframe for '{name}'. Using empty list.")
            kos_dict[name] = []
            continue

        ko_series = df[ko_column_name].tolist()
        parsed_lists = []
        for entry in ko_series:
            parsed = _safe_literal_eval(entry)
            parsed_lists.extend(parsed)

        # unique and sorted
        ko_set_sorted = sorted(set(parsed_lists))
        kos_dict[name] = ko_set_sorted

    return kos_dict


def jaccard(a: set, b: set) -> float:
    """Jaccard similarity for sets. Returns 1.0 when both empty (by convention)."""
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def calculate_similarity_matrix(pattern_kos: Dict[str, List[str]]) -> Tuple[np.ndarray, List[str]]:
    """Return Jaccard similarity matrix (square) and labels.

    Internally converts lists to sets once for speed.
    """
    labels = list(pattern_kos.keys())
    n = len(labels)
    sets = {name: set(kos) for name, kos in pattern_kos.items()}
    matrix = np.zeros((n, n), dtype=float)

    for i, name_i in enumerate(labels):
        si = sets[name_i]
        for j, name_j in enumerate(labels):
            sj = sets[name_j]
            matrix[i, j] = jaccard(si, sj)

    return matrix, labels


def cluster_matrix(matrix: np.ndarray, labels: List[str]) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Cluster using hierarchical average linkage. Returns ordered square matrix, labels, and linkage Z."""
    if matrix.size == 0:
        return matrix, labels, np.array([])

    # distance matrix (0 on diagonal)
    distance_matrix = 1.0 - matrix

    # For clustering, SciPy expects condensed form:
    if len(distance_matrix) == 1:
        # Single item: nothing to cluster
        Z = np.array([])
        leaf_order = np.array([0], dtype=int)
    else:
        condensed = squareform(distance_matrix, checks=True)
        Z = linkage(condensed, method="average")
        leaf_order = leaves_list(Z)

    clustered_matrix = distance_matrix.copy()  # reuse shape
    if clustered_matrix.size:
        clustered_matrix = matrix[np.ix_(leaf_order, leaf_order)]
    clustered_labels = [labels[i] for i in leaf_order]

    return clustered_matrix, clustered_labels, Z


def stat_test(pattern_dict: Dict[str, List[str]], group_labels: Dict[str, str], permutations: int = 5000, seed: int = 5):
    """Run PERMANOVA on the Jaccard distance matrix given valid group labels.

    Validates that groups are provided and that each group has at least 2 samples.
    """
    if not group_labels:
        raise ValueError("No group labels provided for statistical test.")

    matrix, labels = calculate_similarity_matrix(pattern_dict)
    if len(labels) < 2:
        raise ValueError("Need at least 2 samples to run PERMANOVA.")

    # Ensure all labels have group entries
    missing = [lab for lab in labels if lab not in group_labels]
    if missing:
        raise KeyError(f"Missing group labels for samples: {missing}")

    groups = [group_labels[lab] for lab in labels]
    # Check group counts
    counts = pd.Series(groups).value_counts()
    if (counts < 2).any():
        raise ValueError(f"Each group must have at least 2 samples for PERMANOVA. Group sizes:\n{counts.to_dict()}")

    dm = DistanceMatrix(1.0 - matrix, labels)
    return permanova(dm, grouping=groups, permutations=permutations, seed=seed)


def plotting(pattern_dict: Dict[str, List[str]], pattern_name: str, output: str):
    """Create and save heatmap and dendrogram for the given pattern."""
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    matrix, labels = calculate_similarity_matrix(pattern_dict)

    if len(labels) == 0:
        logging.info(f"No samples for pattern '{pattern_name}'; skipping plots.")
        return

    ordered_matrix, ordered_labels, Z = cluster_matrix(matrix=matrix, labels=labels)

    n = len(ordered_labels)
    # Heatmap
    plt.figure(figsize=(max(6, n * 0.5), max(4, n * 0.5)))
    plt.imshow(ordered_matrix, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    plt.xticks(range(n), ordered_labels, rotation=90)
    plt.yticks(range(n), ordered_labels)
    plt.colorbar(label="Jaccard Similarity")
    plt.title("Clustered Jaccard Similarity Heatmap: " + pattern_name)
    plt.tight_layout()
    heatmap_file = output_dir / f"{pattern_name.replace(' ', '')}_GraphComparisons_Heatmap.png"
    plt.savefig(heatmap_file)
    plt.close()
    logging.info(f"Saved heatmap to {heatmap_file}")

    # Dendrogram (only if more than 1 sample)
    if n > 1 and Z.size:
        plt.figure(figsize=(max(6, n * 0.2), 4))
        dendrogram(Z, labels=ordered_labels, leaf_rotation=90)
        plt.title("Hierarchical Clustering Dendrogram: " + pattern_name)
        plt.tight_layout()
        dend_file = output_dir / f"{pattern_name.replace(' ', '')}_GraphComparisons_Dendrogram.png"
        plt.savefig(dend_file)
        plt.close()
        logging.info(f"Saved dendrogram to {dend_file}")
    else:
        logging.info(f"Not enough samples to produce dendrogram for pattern '{pattern_name}' (n={n}).")


def summary(pattern_dict: Dict[str, List[str]], pattern_name: str, stat: bool, groups_dict: Dict[str, str], output: str):
    """Write a summary file including intersection KOs, unique KOs per sample, and optional PERMANOVA results."""
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix, labels = calculate_similarity_matrix(pattern_dict)

    # Intersection across all sets
    all_sets = [set(v) for v in pattern_dict.values()]
    if all_sets:
        intersection_set = set.intersection(*all_sets)
    else:
        intersection_set = set()

    unique_dict = {
        name: sorted(set(ko_list) - intersection_set)
        for name, ko_list in pattern_dict.items()
    }

    summary_file = output_dir / f"{pattern_name.replace(' ', '')}_GraphComparisons_Summary.txt"
    with open(summary_file, "w") as fh:
        fh.write(f"##### SUMMARY FOR PATTERN: {pattern_name} #####\n")
        fh.write(f"Number of compounds shared: {len(intersection_set)}\n")

        if stat:
            try:
                stat_results = stat_test(pattern_dict, groups_dict)
                fh.write("\n##### PERMANOVA RESULTS ######\n")
                fh.write(str(stat_results) + "\n")
            except Exception as e:
                fh.write("\n##### PERMANOVA FAILED ######\n")
                fh.write(f"Reason: {repr(e)}\n")
                logging.warning(f"PERMANOVA failed for pattern '{pattern_name}': {e}")

        for name, unique_list in unique_dict.items():
            fh.write(f"Unique compounds to {name}: {len(unique_list)}\n")

        fh.write("\n\n\n ###### LISTS OF KOs ######\n\n\n")
        fh.write(f"Intersection KOs: {sorted(intersection_set)}\n\n")

        for name, unique_list in unique_dict.items():
            fh.write(f"Unique KOs for {name}: {unique_list}\n\n")

    logging.info(f"Saved summary to {summary_file}")


def main():
    parser = argparse.ArgumentParser(description="Compare graph results across samples using KOs and Jaccard similarity.")
    parser.add_argument("-m", "--metadata", required=True, help="Metadata CSV containing file paths, names, and groups")
    parser.add_argument("-p", "--paths", required=True, help="Name of column containing file paths")
    parser.add_argument("-n", "--names", required=True, help="Name of column containing names of graphs (e.g., sampleID)")
    parser.add_argument("-s", "--stat_test", action="store_true", help="If statistical test for group comparison wanted include this parameter")
    parser.add_argument("-g", "--groups", help="Name of column containing group labels", default="")
    parser.add_argument("-o", "--output", required=True, help="Output directory for plots and summary files")
    parser.add_argument("--ko_column", help="Name of KOs column in graph CSVs (default: 'KOs')", default="KOs")
    args = parser.parse_args()

    if args.stat_test and not args.groups:
        parser.error("Statistical test requested but no groups column provided (-g/--groups).")

    paths, names, group_dict = csv_to_inputs(metadata=args.metadata, paths_col=args.paths, names_col=args.names, groups_col=args.groups if args.stat_test else "")

    graphs_dict = get_graphs(paths=paths, names=names)
    food_microbe_dict, food_both_dict, both_both_dict = subset_graphs(graph_dict=graphs_dict)

    food_microbe_kos = get_kos(food_microbe_dict, ko_column_name=args.ko_column)
    food_both_kos = get_kos(food_both_dict, ko_column_name=args.ko_column)
    both_both_kos = get_kos(both_both_dict, ko_column_name=args.ko_column)

    patterns = [food_microbe_kos, food_both_kos, both_both_kos]
    pattern_names = ["Food to Microbe", "Food to Both", "Both to Both"]

    for pat_dict, pat_name in zip(patterns, pattern_names):
        plotting(pat_dict, pat_name, args.output)
        summary(pattern_dict=pat_dict, pattern_name=pat_name, stat=args.stat_test, groups_dict=group_dict, output=args.output)


if __name__ == "__main__":
    main()