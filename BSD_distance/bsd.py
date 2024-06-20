'''
This script calculates the Branch Score Distance (BSD) between two phylogenetic trees provided in Newick format. 
Depending on the overlap of the leaf sets of the two trees, it calculates either the BSD or the BSD(-) value.
'''

import argparse
from ete3 import Tree
from itertools import combinations
import math

# Step 1: Parsing and Initial Analysis
def parse_newick(newick_str):
    return Tree(newick_str, format=1)

def find_common_leaves(tree1, tree2):
    leaves1 = set(leaf.name for leaf in tree1.get_leaves())
    leaves2 = set(leaf.name for leaf in tree2.get_leaves())
    return leaves1.intersection(leaves2)

# Step 2: Distance Calculations
def calculate_pairwise_distances(tree, leaves):
    distances = {}
    for leaf1, leaf2 in combinations(leaves, 2):
        distance = tree.get_distance(leaf1, leaf2)
        distances[(leaf1, leaf2)] = distance
    return distances

# Step 3: Calculate the Branch Score Distance (BSD)
def calculate_BSD(tree1, tree2, leaves):
    def squared_distance_sum(t1, t2, leaves):
        sum_sq_distance = 0
        for leaf1, leaf2 in combinations(leaves, 2):
            d1 = t1.get_distance(leaf1, leaf2)
            d2 = t2.get_distance(leaf1, leaf2)
            sum_sq_distance += (d1 - d2) ** 2
        return sum_sq_distance
    return math.sqrt(squared_distance_sum(tree1, tree2, leaves))

# Function to prune a tree to only contain common leaves
def prune_to_common_leaves(tree, common_leaves):
    pruned_tree = tree.copy()
    pruned_tree.prune(common_leaves)
    return pruned_tree

def main(input_file, output_file):
    with open(input_file, 'r') as file:
        trees = [parse_newick(line.strip()) for line in file if line.strip()]

    if len(trees) != 2:
        print("Please provide exactly two trees in the input file.")
        return

    tree1 = trees[0]
    tree2 = trees[1]
    common_leaves = find_common_leaves(tree1, tree2)

    if not common_leaves:
        print("The BSD distance between input trees cannot be computed because these trees have no common leaves.")
        return

    leaves1 = set(leaf.name for leaf in tree1.get_leaves())
    leaves2 = set(leaf.name for leaf in tree2.get_leaves())

    with open(output_file, 'w') as out:
        if leaves1 == leaves2:
            # Calculate BSD
            bsd = calculate_BSD(tree1, tree2, common_leaves)
            bsd_rounded = round(bsd, 4)
            out.write(f"BSD(+): {bsd_rounded}\n")
        else:
            # Prune trees and calculate BSD-minus
            pruned_tree1 = prune_to_common_leaves(tree1, common_leaves)
            pruned_tree2 = prune_to_common_leaves(tree2, common_leaves)
            bsd_minus = calculate_BSD(pruned_tree1, pruned_tree2, common_leaves)
            bsd_minus_rounded = round(bsd_minus, 4)
            out.write(f"BSD(-): {bsd_minus_rounded}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate BSD distances between two phylogenetic trees.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file containing two Newick trees')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file to write results')
    
    args = parser.parse_args()
    main(args.input, args.output)
