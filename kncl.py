# Please note that this version is not final and is under development

import argparse
from ete3 import Tree
import math
from itertools import combinations

# Helper functions
def precompute_descendants(node, distinct_leaves):
    if node.is_leaf():
        node.add_feature("descendants_distinct", node.name in distinct_leaves)
    else:
        all_distinct = True
        for child in node.children:
            if not child.descendants_distinct:
                all_distinct = False
                break
        node.add_feature("descendants_distinct", all_distinct)

def findSD(tree, distinct_leaves):
    for node in tree.traverse("postorder"):
        precompute_descendants(node, distinct_leaves)

    subtree_roots = set()
    visited_leaves = set()

    for leaf_name in distinct_leaves:
        if leaf_name in visited_leaves:
            continue

        leaf_node = tree & leaf_name
        current_node = leaf_node
        subtree_root = None

        while current_node:
            if current_node.descendants_distinct:
                subtree_root = current_node
            current_node = current_node.up

        if subtree_root and subtree_root not in subtree_roots:
            subtree_roots.add(subtree_root)
            visited_leaves.update(get_leaves(subtree_root))

    return subtree_roots

def get_leaves(node):
    return set(leaf.name for leaf in node)

def get_subtree_newick_with_branch_lengths(node):
    return node.write(format=1)

def InsertTempLeaves(tree, target_leaf, new_leaf_base_name, new_length, dist, inserted_leaves, tolerance=1e-10):
    # Operate directly on 'tree'
    target_node = tree.search_nodes(name=target_leaf)[0]
    insertion_points = []
    visited_nodes = set()

    # Label internal nodes and branches for easier tracking
    internal_node_counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and not node.name:
            node.name = f"Node{internal_node_counter}"
            internal_node_counter += 1

    def robust_insert_leaf_at_node(current_node, insert_distance, previous_node, original_branch_distance, toward_root=False):
        excess_length = original_branch_distance - insert_distance

        if excess_length < 0:
            excess_length = 0

        # Handle traversal toward the root by ensuring correct branch selection
        if toward_root:
            temp = current_node
            current_node = previous_node
            previous_node = temp

        # Detach the previous node from its parent (the node leading to the root)
        parent = previous_node.up
        if parent:
            previous_node.detach()

        if parent is None:
            new_internal_node = tree.add_child(dist=excess_length)
            current_node.detach()
            new_internal_node.add_child(current_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_leaf_name)
            visited_nodes.add(new_internal_node)
        else:
            new_internal_node = parent.add_child(dist=excess_length)
            new_internal_node.add_child(previous_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_leaf_name)
            visited_nodes.add(new_internal_node)

        return True

    def insert_leaf_at_terminal(current_node, insert_distance):
        if current_node.name in inserted_leaves:
            return False
        excess_length = current_node.dist - insert_distance
        if excess_length < 0:
            excess_length = 0

        parent = current_node.up
        if parent:
            current_node.detach()
            new_internal_node = parent.add_child(dist=excess_length)
            new_internal_node.add_child(current_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_leaf_name)
            visited_nodes.add(new_internal_node)
        else:
            return False

        return True

    def bfs(node, accumulated_distance):
        queue = [(node, accumulated_distance, None, 0, [], False)]
        while queue:
            current_node, current_dist, prev_node, prev_dist, path, toward_root = queue.pop(0)
            if current_node in visited_nodes or 'temp' in current_node.name or current_node.name in inserted_leaves:
                continue
            visited_nodes.add(current_node)
            current_path = path + [current_node.name]

            if round(current_dist, 8) >= dist:
                insert_distance = round(current_dist, 8) - round(dist, 8)
                if abs(insert_distance) < tolerance:
                    insert_distance = 0
                if insert_distance == 0:
                    if not robust_insert_leaf_at_node(current_node, insert_distance, prev_node, current_node.dist, toward_root):
                        return
                elif current_node.is_leaf():
                    if not insert_leaf_at_terminal(current_node, insert_distance):
                        return
                else:
                    if not robust_insert_leaf_at_node(prev_node, prev_dist - insert_distance, current_node, prev_dist, toward_root):
                        return
                continue

            for child in current_node.children:
                if child not in visited_nodes and child.name not in inserted_leaves:
                    queue.append((child, current_dist + child.dist, current_node, child.dist, current_path, False))

            if current_node.up and current_node.up not in visited_nodes and current_node.up.name not in inserted_leaves:
                queue.append((current_node.up, current_dist + current_node.dist, current_node, current_node.dist, current_path, True))

    if dist <= target_node.dist:
        insert_leaf_at_terminal(target_node, dist)
    else:
        bfs(target_node, 0)

    return insertion_points  # Return names instead of node objects

def find_farthest_leaf(tree, start, temporary_leaves):
    max_distance = 0
    farthest_leaf = start
    for leaf_name in temporary_leaves:
        if leaf_name != start.name:
            leaf = tree & leaf_name
            distance = tree.get_distance(start, leaf)
            if distance > max_distance:
                max_distance = distance
                farthest_leaf = leaf
    return farthest_leaf, max_distance

def find_path(leaf1, leaf2):
    path1 = []
    node = leaf1
    while node:
        path1.append(node)
        node = node.up

    path2 = []
    node = leaf2
    while node:
        path2.append(node)
        node = node.up

    lca = None
    for n1 in path1:
        if n1 in path2:
            lca = n1
            break

    path = []
    branch_lengths = []
    for node in path1:
        path.append(node)
        if node == lca:
            break

    lca_index = path2.index(lca)
    reversed_path2 = list(reversed(path2[:lca_index]))
    path.extend(reversed_path2)

    for i in range(1, len(path)):
        branch_lengths.append(path[i - 1].get_distance(path[i]))

    return path, branch_lengths

def compute_midpoint(tree, temporary_leaves):
    start_name = next(iter(temporary_leaves))
    start = tree & start_name
    leaf1, dist1 = find_farthest_leaf(tree, start, temporary_leaves)
    leaf2, dist2 = find_farthest_leaf(tree, leaf1, temporary_leaves)
    path, branch_lengths = find_path(leaf1, leaf2)
    total_distance = dist2
    half_distance = round(total_distance / 2, 8)

    cumulative_distance = 0
    prev_node = None
    for i, node in enumerate(path):
        if i > 0:
            cumulative_distance += round(branch_lengths[i - 1], 8)
        if cumulative_distance >= half_distance:
            excess = round(cumulative_distance - half_distance, 8)
            prev_node = path[i - 1]
            return prev_node, node, excess, half_distance, branch_lengths[i - 1]
        elif cumulative_distance == half_distance:
            prev_node = path[i - 1]
            return prev_node, node, 0, half_distance, branch_lengths[i - 1]

def insert_midpoint_and_new_subtree(tree, prev_node, curr_node, excess, subtree, branch_length, original_dist):
    if excess == 0:
        # Attach the subtree directly to curr_node
        new_subtree = subtree.copy()
        new_subtree.dist = branch_length
        curr_node.add_child(new_subtree)
        return tree

    if curr_node in prev_node.get_ancestors():
        distance_to_midpoint = round(excess, 8)
        distance_from_midpoint_to_leaf = round(original_dist - excess, 8)
    else:
        distance_to_midpoint = round(original_dist - excess, 8)
        distance_from_midpoint_to_leaf = round(excess, 8)

    if distance_to_midpoint < 0 or distance_from_midpoint_to_leaf < 0:
        raise ValueError("Negative distance encountered. Check the calculation logic.")

    new_node = Tree()
    new_node.dist = distance_to_midpoint

    if curr_node in prev_node.get_ancestors():
        parent = prev_node.up
        child = prev_node
    else:
        parent = prev_node
        child = curr_node

    parent.remove_child(child)
    parent.add_child(new_node)
    new_node.add_child(child)
    child.dist = distance_from_midpoint_to_leaf

    # Now add the subtree
    new_subtree = subtree.copy()
    new_subtree.dist = branch_length
    new_node.add_child(new_subtree)

    return tree

def remove_temporary_leaves(tree, temporary_leaves):
    def collapse_single_child_nodes(node):
        while not node.is_leaf() and len(node.children) == 1:
            child = node.children[0]
            child.dist += node.dist
            parent = node.up
            if parent:
                parent.remove_child(node)
                parent.add_child(child)
                node = parent
            else:
                # Node is root
                child.up = None
                tree = child  # Update tree reference
                node = child
    for leaf_name in temporary_leaves:
        leaf_nodes = tree.search_nodes(name=leaf_name)
        if leaf_nodes:
            leaf = leaf_nodes[0]
            parent = leaf.up
            if parent:
                parent.remove_child(leaf)
                collapse_single_child_nodes(parent)
            else:
                # Leaf is root
                tree = None
    return tree  # Return the possibly updated tree

def clear_internal_node_names(tree):
    for node in tree.traverse():
        if not node.is_leaf():
            node.name = ''

def kNCL(T1, T2, k):
    CL = set(T1.get_leaf_names()) & set(T2.get_leaf_names())
    if len(CL) < 3:
        raise ValueError("The input trees must have at least 3 common leaves.")
    if k < 2 or k > len(CL):
        raise ValueError("The value of k must be between 2 and the number of common leaves.")

    def adjust_rate(T1, T2):
        sum_T1 = sum(T1.get_distance(l1, l2) for i, l1 in enumerate(CL) for l2 in list(CL)[i + 1:])
        sum_T2 = sum(T2.get_distance(l1, l2) for i, l1 in enumerate(CL) for l2 in list(CL)[i + 1:])
        return sum_T1 / sum_T2 if sum_T2 else 1

    r12 = adjust_rate(T1, T2)  # Adjusting from T2 to T1
    r21 = adjust_rate(T2, T1)  # Adjusting from T1 to T2

    SD1 = findSD(T1, set(T1.get_leaf_names()) - CL)
    SD2 = findSD(T2, set(T2.get_leaf_names()) - CL)

    inserted_leaves = set()  # Track inserted leaves to ignore them in future iterations

    def process_tree(target_tree, source_tree, subtrees_to_insert, rate, k):
        for a in subtrees_to_insert:
            if not a.name:
                a.name = "subtree_" + str(len(subtrees_to_insert))  # Assign a name to unnamed subtrees

            # Make a copy of the subtree to avoid modifying the original
            adjusted_subtree = a.copy()

            # Adjust all branch lengths in the copied subtree
            for node in adjusted_subtree.traverse():
                node.dist *= rate

            NCL = sorted(CL, key=lambda l: source_tree.get_distance(a, source_tree & l))[:k]
            TL = set()
            for lc in NCL:
                if target_tree.search_nodes(name=lc):
                    lc_node = target_tree & lc
                else:
                    raise ValueError(f"Common leaf '{lc}' not found in target_tree.")
                lc_node_source = source_tree & lc
                rc = sum(target_tree.get_distance(lc_node, l) for l in CL) / sum(source_tree.get_distance(lc_node_source, l) for l in CL)
                dp = (source_tree.get_distance(a, lc_node_source) - a.dist) * rc
                temp_leaves = InsertTempLeaves(target_tree, lc, "temp", adjusted_subtree.dist, dp, inserted_leaves)
                TL.update(temp_leaves)
                inserted_leaves.update(temp_leaves)  # Add the inserted leaves to the set


            if not TL:
                continue

            if len(TL) == 1:
                single_leaf_name = next(iter(TL))
                single_leaf = target_tree & single_leaf_name
                parent = single_leaf.up
                branch_length = single_leaf.dist

                # Remove the temporary leaf
                parent.remove_child(single_leaf)

                # Attach the adjusted subtree to the parent
                adjusted_subtree.dist = branch_length
                parent.add_child(adjusted_subtree)

            else:
                try:
                    prev_node, curr_node, excess, _, original_dist = compute_midpoint(target_tree, TL)
                    target_tree = insert_midpoint_and_new_subtree(target_tree, prev_node, curr_node, excess, adjusted_subtree, adjusted_subtree.dist, original_dist)
                except Exception as e:
                    print("Error encountered during midpoint insertion:")
                    print(f"Nodes involved: {TL}")
                    print(f"Error: {str(e)}")

            target_tree = remove_temporary_leaves(target_tree, TL)
        return target_tree

    # Process trees with the corrected adjustment rates
    T1_completed = process_tree(T1, T2, SD2, r12, k)

    T2_completed = process_tree(T2, T1, SD1, r21, k)

    clear_internal_node_names(T1_completed)
    clear_internal_node_names(T2_completed)

    return T1_completed, T2_completed

def parse_input_file(input_file):
    with open(input_file, 'r') as file:
        tree_lines = file.readlines()
    trees = [Tree(tree_str.strip(), format=1) for tree_str in tree_lines]
    return trees

def write_output_file(output_file, results):
    with open(output_file, 'w') as file:
        for result in results:
            file.write(result + '\n')

def squared_distance_sum(t1, t2, leaves):
    sum_sq_distance = 0
    for leaf1, leaf2 in combinations(leaves, 2):
        d1 = t1.get_distance(leaf1, leaf2)
        d2 = t2.get_distance(leaf1, leaf2)
        sum_sq_distance += (d1 - d2) ** 2
    return sum_sq_distance

def BSD(T1, T2, k):
    # Get the leaves from the original input trees (before completion)
    leaves1 = set(leaf.name for leaf in T1.get_leaves())
    leaves2 = set(leaf.name for leaf in T2.get_leaves())

    # Find the common leaves between the two original trees
    common_leaves = leaves1.intersection(leaves2)
    if len(common_leaves) < 3:
        return None, None, None, None

    # Complete both trees using the k-NCL function
    T1_completed, T2_completed = kNCL(T1, T2, k)

    # Get the leaves of the completed trees
    leaves_completed = set(leaf.name for leaf in T1_completed.get_leaves())

    # Calculate BSD(+) over the completed trees and the leafset of T1_completed
    bsd_plus = math.sqrt(squared_distance_sum(T1_completed, T2_completed, leaves_completed))

    # Calculate BSD(-) over the common leaves of the original trees
    bsd_minus = math.sqrt(squared_distance_sum(T1, T2, common_leaves))

    # Return BSD distances and the completed trees in Newick format
    return bsd_plus, bsd_minus, T1_completed.write(format=1), T2_completed.write(format=1)

def main():
    parser = argparse.ArgumentParser(description="Run k-NCL algorithm on Newick trees.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input file with trees in Newick format")
    parser.add_argument('k', type=int, help="Integer value of k for k-NCL algorithm")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file to save results")

    args = parser.parse_args()

    # Parse the input trees
    trees = parse_input_file(args.input)
    if len(trees) < 2:
        raise ValueError("The input file must contain at least two trees.")

    # Run the k-NCL algorithm and calculate BSD for each pair of trees
    results = []
    for i in range(len(trees)):
        for j in range(i + 1, len(trees)):
            T1 = trees[i]
            T2 = trees[j]
            bsd_plus, bsd_minus, T1_completed_newick, T2_completed_newick = BSD(T1, T2, args.k)
            if bsd_plus is not None and bsd_minus is not None:
                result = (f"Tree pair {i + 1} and {j + 1}:\n"
                          f"BSD(+) = {bsd_plus:.4f}, BSD(-) = {bsd_minus:.4f}\n"
                          f"Completed Tree 1:\n{T1_completed_newick}\n"
                          f"Completed Tree 2:\n{T2_completed_newick}")
                results.append(result)
            else:
                results.append(f"Tree pair {i + 1} and {j + 1}: Tree completion cannot be performed on these trees. Check their common leaves.")

    # Write the output to a file
    write_output_file(args.output, results)

if __name__ == "__main__":
    main()
