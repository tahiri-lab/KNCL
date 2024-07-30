# Finding the midpoint among a set of temporary leaves TL and inserting a new leaf at the midpoint

from ete3 import Tree

def find_farthest_leaf(tree, start):
    # Retrieve all leaves in the tree
    leaves = tree.get_leaves()
    max_distance = 0
    farthest_leaf = start

    for leaf in leaves:
        if leaf is not start:  # Ensure we don't measure distance to itself
            distance = tree.get_distance(start, leaf)
            if distance > max_distance:
                max_distance = distance
                farthest_leaf = leaf

    print(f"Farthest leaf from {start.name}: {farthest_leaf.name} at distance {max_distance}")
    return farthest_leaf, max_distance

def find_path(tree, leaf1, leaf2):
    path = []
    # Get the path from each leaf to the root and find the common ancestor
    node = leaf1
    path1 = []
    while node:
        path1.append(node)
        node = node.up
    node = leaf2
    path2 = []
    while node:
        path2.append(node)
        node = node.up

    lca = set(path1) & set(path2)
    for node in path1:
        if node in lca:
            break
        path.append(node)
    path.append(node)  # Append the lowest common ancestor
    for node in reversed(path2):
        if node in lca:
            continue
        path.append(node)

    path_names = [n.name for n in path]
    print(f"Diameter path: {' -> '.join(path_names)}")
    return path

def compute_midpoint(tree, temporary_leaves):
    start = next(iter(temporary_leaves))
    leaf1, dist1 = find_farthest_leaf(tree, start)
    leaf2, dist2 = find_farthest_leaf(tree, leaf1)
    path = find_path(tree, leaf1, leaf2)
    total_distance = dist2
    half_distance = total_distance / 2

    cumulative_distance = 0
    for node in path:
        cumulative_distance += node.dist
        if cumulative_distance > half_distance:
            excess = cumulative_distance - half_distance
            return node, path, node.dist - excess

def insert_midpoint_and_new_leaf(tree, midpoint_node, excess, new_leaf_name, branch_length):
    # Calculate distances
    parent = midpoint_node.up
    original_dist = midpoint_node.dist
    distance_to_midpoint = original_dist - excess
    distance_from_midpoint_to_leaf = excess

    # Ensure shorter distance is towards the midpoint
    if distance_to_midpoint > distance_from_midpoint_to_leaf:
        distance_to_midpoint, distance_from_midpoint_to_leaf = distance_from_midpoint_to_leaf, distance_to_midpoint

    # Insert new midpoint node
    new_node = Tree(name="midpoint")
    new_node.dist = distance_to_midpoint

    # Reconfigure the parent and child relationships
    parent.remove_child(midpoint_node)
    parent.add_child(new_node)
    midpoint_node.dist = distance_from_midpoint_to_leaf
    new_node.add_child(midpoint_node)

    # Add the new leaf
    new_leaf = Tree(name=new_leaf_name)
    new_leaf.dist = branch_length
    new_node.add_child(new_leaf)

    print(f"Inserted '{new_leaf_name}' at midpoint with excess distance: {excess}")
    print(f"Branch lengths: {distance_to_midpoint} to midpoint, {distance_from_midpoint_to_leaf} to leaf")

    return tree

# Example Usage
newick = "((A:1.5,B:1.2)Node1:1.5,(C:1.9,D:0.1)Node2:0.5);"
tree = Tree(newick, format=1)
print("Original tree:")
print(tree)

temporary_leaves = {tree & "A", tree & "C"}
new_leaf_name = "New_leaf"
branch_length = 1.5

midpoint_node, path, excess = compute_midpoint(tree, temporary_leaves)
tree = insert_midpoint_and_new_leaf(tree, midpoint_node, excess, new_leaf_name, branch_length)

print("Updated tree:")
print(tree)
print(tree.write(format=1))
