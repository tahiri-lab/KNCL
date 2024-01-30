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

def find_distinct_leaves(tree, common_leaves):
    return set(leaf.name for leaf in tree.get_leaves()) - common_leaves

def find_maximal_distinct_leaf_subtrees(tree, distinct_leaves):
    """Find maximal distinct-leaf subtrees in a tree."""
    MDS = []
    subtree_roots = set()

    for leaf_name in distinct_leaves:
        leaf = tree & leaf_name  # Find the leaf node
        current_node = leaf

        # Traverse up the tree until a node is reached that is not part of the maximal subtree
        while current_node:
            leaves_under_current_node = set(node.name for node in current_node.get_leaves())
            if leaves_under_current_node.issubset(distinct_leaves):
                subtree_root = current_node
                current_node = current_node.up
            else:
                break

        # Add the subtree root if it's not already included
        if subtree_root not in subtree_roots:
            subtree_roots.add(subtree_root)
            MDS.append(set(node.name for node in subtree_root.get_leaves()))

    return MDS

# Step 2: Distance Calculations and Adjustment Rates
def calculate_pairwise_distances(tree, leaves):
    distances = {}
    for leaf1, leaf2 in combinations(leaves, 2):
        distance = tree.get_distance(leaf1, leaf2)
        distances[(leaf1, leaf2)] = distance
    return distances

def calculate_leaf_rates(distances1, distances2, common_leaves):
    rates = {}
    for leaf in common_leaves:
        sum_distances1 = sum(distances1.get((leaf, other_leaf), 0) for other_leaf in common_leaves if other_leaf != leaf)
        sum_distances2 = sum(distances2.get((leaf, other_leaf), 0) for other_leaf in common_leaves if other_leaf != leaf)
        rates[leaf] = sum_distances1 / sum_distances2 if sum_distances2 != 0 else 1
    return rates

# Step 4: Adding Elements to Trees
def calculate_cutback_distances(source_tree, element, common_leaves, k):
    if element not in [n.name for n in source_tree.get_leaves()]:
        #print(f"Element '{element}' not found in the source tree.")
        return {}
    
    element_node = source_tree & element
    distances = {leaf: source_tree.get_distance(leaf, element_node) - element_node.dist 
                 for leaf in common_leaves if leaf in [n.name for n in source_tree.get_leaves()]}
    sorted_distances = dict(sorted(distances.items(), key=lambda item: item[1]))
    return dict(list(sorted_distances.items())[:k])

def adjust_distances_by_rates(distances, rates):
    return {leaf: distances[leaf] * rates.get(leaf, 1) for leaf in distances}

def find_longest_path(tree, nodes):
    max_distance = 0
    node_pair = (None, None)
    for node1, node2 in combinations(nodes, 2):
        distance = tree.get_distance(node1, node2)
        if distance > max_distance:
            max_distance = distance
            node_pair = (node1, node2)
    return node_pair

def get_path_between_nodes(node1, node2):
    if not node1 or not node2:
        # If either node is None, return an empty path
        return []

    path_to_root = set()
    current = node1
    while current:
        path_to_root.add(current)
        current = current.up

    path = []
    current = node2
    while current and current not in path_to_root:
        path.append(current)
        current = current.up

    if current is None:
        # If there's no intersection, the nodes are not part of the same tree
        return []

    while node1 and node1 != current:
        path.append(node1)
        node1 = node1.up

    return path

def find_midpoint_on_path(path):
    if not path:  # Check if the path is empty
        #print("No path found for midpoint calculation.")
        return None

    mid_index = len(path) // 2
    return path[mid_index]

def add_new_element_to_tree(tree, element, adjusted_distances):
    #print(f"Adding new element '{element}' to the tree")

    for common_leaf, distance in adjusted_distances.items():
        if common_leaf in [n.name for n in tree.get_leaves()]:
            common_leaf_node = tree & common_leaf

            # Clone the subtree rooted at the common leaf node
            cloned_subtree = common_leaf_node.copy()

            # Create new internal node and add the common leaf and new element as children
            new_internal_node = Tree()
            new_internal_node.add_child(Tree(name=element, dist=distance))
            new_internal_node.add_child(cloned_subtree)

            # Attach the new internal node to the parent of the common leaf node
            if common_leaf_node.up:
                parent = common_leaf_node.up
                parent.add_child(new_internal_node)
                parent.remove_child(common_leaf_node)
            else:
                # If common_leaf_node is the root, make new_internal_node the new root
                tree.set_outgroup(new_internal_node)

            #print(f"New element '{element}' added near '{common_leaf}'")
            break

def add_element_to_tree(target_tree, source_tree, element, common_leaves, leaf_rates, k):
    cutback_distances = calculate_cutback_distances(source_tree, element, common_leaves, k)
    adjusted_distances = adjust_distances_by_rates(cutback_distances, leaf_rates)

    if adjusted_distances:
        add_new_element_to_tree(target_tree, element, adjusted_distances)
    else:
        print(f"No adjusted distances for '{element}', skipping addition.")

# Step 5: Calculate the Branch Score Distance (BSD)
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

def join_trees_with_new_root(tree1, tree2):
    """Joins two trees with a new root node if they have no common leaves."""
    new_root = Tree()  # Create a new empty root node
    new_root.add_child(tree1)  # Add the first tree as a child
    new_root.add_child(tree2)  # Add the second tree as a child
    return new_root

def main(input_file, k, output_file):
    with open(input_file, 'r') as file:
        trees = [parse_newick(line.strip()) for line in file if line.strip()]

    with open(output_file, 'w') as out:
        tree1 = trees[0]
        tree1_copy = tree1.copy()
        for tree2 in trees[1:]:
            tree2_copy = tree2.copy()
            common_leaves = find_common_leaves(tree1_copy, tree2_copy)
            if not common_leaves:
                completed_tree = join_trees_with_new_root(tree1_copy, tree2_copy)
                completed_tree_str = completed_tree.write(format=1)
                out.write(f"Completed Tree (joined with a new root):\n{completed_tree_str}\n")
                bsd = bsd_minus = 'N/A'  # BSD is not applicable in this case
            else:
                # Proceed with normal completion process
                distances_tree1 = calculate_pairwise_distances(tree1_copy, common_leaves)
                distances_tree2 = calculate_pairwise_distances(tree2_copy, common_leaves)
                leaf_rates = calculate_leaf_rates(distances_tree1, distances_tree2, common_leaves)

                # Perform tree completion (Step 4)
                distinct_leaves_tree1 = find_distinct_leaves(tree1_copy, common_leaves)
                for leaf in distinct_leaves_tree1:
                    add_element_to_tree(tree2_copy, tree1, leaf, common_leaves, leaf_rates, k)

                # Repeat for distinct leaves from tree2 to tree1
                distinct_leaves_tree2 = find_distinct_leaves(tree2_copy, common_leaves)
                for leaf in distinct_leaves_tree2:
                    add_element_to_tree(tree1_copy, tree2, leaf, common_leaves, leaf_rates, k)

                completed_tree1_str = tree1_copy.write(format=1)
                completed_tree2_str = tree2_copy.write(format=1)
                # Calculate BSD
                bsd = calculate_BSD(tree1_copy, tree2_copy, common_leaves)
                # Prune trees and calculate BSD-minus
                pruned_tree1 = prune_to_common_leaves(tree1_copy, common_leaves)
                pruned_tree2 = prune_to_common_leaves(tree2_copy, common_leaves)
                pruned_tree1_str = pruned_tree1.write(format=1)
                pruned_tree2_str = pruned_tree2.write(format=1)
                bsd_minus = calculate_BSD(pruned_tree1, pruned_tree2, common_leaves)

            # Write the results to the output file
            out.write(f"Completed Tree 1:\n{completed_tree1_str}\n")
            out.write(f"Completed Tree 2:\n{completed_tree2_str}\n")
            out.write(f"BSD(+): {bsd}\n")
            out.write(f"Pruned Tree 1:\n{pruned_tree1_str}\n")
            out.write(f"Pruned Tree 2:\n{pruned_tree2_str}\n")
            out.write(f"BSD(-): {bsd_minus}\n\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Complete phylogenetic trees and calculate BSD distances between them.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file containing Newick trees')
    parser.add_argument('k', type=int, help='Value of k for the k-NCL algorithm')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file to write results')
    
    args = parser.parse_args()
    main(args.input, args.k, args.output)