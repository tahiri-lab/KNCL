# Finding the midpoint among a set of temporary leaves TL and inserting a new leaf at the midpoint

from ete3 import Tree

def find_farthest_leaf(tree, start, temporary_leaves):
    distances = {}
    farthest_leaf = start
    max_distance = 0

    for leaf in temporary_leaves:
        distance = tree.get_distance(start, leaf)
        distances[leaf] = distance
        if distance > max_distance:
            max_distance = distance
            farthest_leaf = leaf

    print(f"Farthest leaf from {start.name}: {farthest_leaf.name}")
    return farthest_leaf, distances

def find_path(leaf1, leaf2):
    common_ancestor = leaf1.get_common_ancestor(leaf2)
    path1 = []
    current = leaf1
    while current != common_ancestor:
        path1.append(current)
        current = current.up
    path1.append(common_ancestor)
    path2 = []
    current = leaf2
    while current != common_ancestor:
        path2.append(current)
        current = current.up
    path2.reverse()
    path = path1 + path2[1:]
    path_names = [leaf.name for leaf in path]
    print(f"Diameter path: {' -> '.join(path_names)}")
    return path

def compute_midpoint(tree, temporary_leaves):
    start = next(iter(temporary_leaves))
    print(f"Starting leaf: {start.name}")
    leaf1, _ = find_farthest_leaf(tree, start, temporary_leaves)
    leaf2, distances = find_farthest_leaf(tree, leaf1, temporary_leaves)
    diameter_path = find_path(leaf1, leaf2)
    total_distance = distances[leaf2]
    half_distance = total_distance / 2

    cumulative_distance = 0
    for i, node in enumerate(diameter_path):
        if i > 0:
            cumulative_distance += diameter_path[i-1].dist
        if cumulative_distance >= half_distance:
            break

    midpoint_position = diameter_path[i]
    print(f"Midpoint position: {midpoint_position.name}")
    return midpoint_position, diameter_path, half_distance

def insert_midpoint_and_new_leaf(tree, diameter_path, half_distance, new_leaf_name, branch_length):
    cumulative_distance = 0
    for i, node in enumerate(diameter_path):
        cumulative_distance += node.dist
        if cumulative_distance >= half_distance:
            break

    midpoint_position = diameter_path[i]
    parent = midpoint_position.up

    # Calculate the exact distance to the midpoint for accurate split
    distance_to_midpoint = cumulative_distance - half_distance
    print(f"Distance to midpoint: {distance_to_midpoint}")

    new_node = Tree(name="midpoint")
    new_node.dist = distance_to_midpoint
    midpoint_position.dist -= distance_to_midpoint

    parent.add_child(new_node)
    new_node.add_child(midpoint_position.detach())

    print(f"Inserting new leaf '{new_leaf_name}' with branch length {branch_length} at midpoint")
    new_leaf = Tree(name=new_leaf_name)
    new_leaf.dist = branch_length
    new_node.add_child(new_leaf)

    return tree

# Example

newick = "((A:1.5,B:1.2):1.5,(C:1.9,D:0.1):0.5);"
tree = Tree(newick, format=1)

print("Original tree:")
print(tree)

temporary_leaves = {tree & "A", tree & "C", tree & "D"}

new_leaf_name = "New_leaf"
branch_length = 1.5

midpoint_position, diameter_path, half_distance = compute_midpoint(tree, temporary_leaves)
tree = insert_midpoint_and_new_leaf(tree, diameter_path, half_distance, new_leaf_name, branch_length)

print("Updated tree:")
print(tree)
print(tree.write(format=1))

# Updated approach with debugging statements

from ete3 import Tree

def label_internal_nodes(tree):
    for i, node in enumerate(tree.traverse("postorder")):
        if not node.is_leaf() and not node.name:
            node.name = f"Node{i}"

def find_farthest_leaf(tree, start, temporary_leaves):
    distances = {}
    farthest_leaf = start
    max_distance = 0

    for leaf in temporary_leaves:
        distance = tree.get_distance(start, leaf)
        distances[leaf] = distance
        if distance > max_distance:
            max_distance = distance
            farthest_leaf = leaf

    print(f"Farthest leaf from {start.name}: {farthest_leaf.name}, Distance: {max_distance}")
    return farthest_leaf, distances

def find_path(leaf1, leaf2):
    common_ancestor = leaf1.get_common_ancestor(leaf2)
    path1 = []
    current = leaf1
    while current != common_ancestor:
        path1.append(current)
        current = current.up
    path1.append(common_ancestor)
    path2 = []
    current = leaf2
    while current != common_ancestor:
        path2.append(current)
        current = current.up
    path2.reverse()
    path = path1 + path2[1:]
    path_names = [leaf.name for leaf in path]
    print(f"Diameter path: {' -> '.join(path_names)}")
    return path

def compute_midpoint(tree, temporary_leaves):
    start = next(iter(temporary_leaves))
    print(f"Starting leaf: {start.name}")
    leaf1, _ = find_farthest_leaf(tree, start, temporary_leaves)
    leaf2, distances = find_farthest_leaf(tree, leaf1, temporary_leaves)
    diameter_path = find_path(leaf1, leaf2)
    total_distance = distances[leaf2]
    half_distance = total_distance / 2
    print(f"Total distance between {leaf1.name} and {leaf2.name}: {total_distance}, Half distance: {half_distance}")

    cumulative_distance = 0
    for i, node in enumerate(diameter_path):
        if i > 0:
            cumulative_distance += diameter_path[i-1].dist
        print(f"Traversing Node: {node.name}, Dist: {node.dist}, Cumulative distance: {cumulative_distance}")
        if cumulative_distance + node.dist >= half_distance:
            break

    if cumulative_distance + node.dist > half_distance:
        distance_to_midpoint = cumulative_distance + node.dist - half_distance
        midpoint_position = diameter_path[i]
    else:
        distance_to_midpoint = half_distance - cumulative_distance
        midpoint_position = diameter_path[i + 1]

    print(f"Midpoint position: {midpoint_position.name}, Cumulative distance: {cumulative_distance}, Distance to midpoint: {distance_to_midpoint}")
    return midpoint_position, diameter_path, half_distance, distance_to_midpoint

def insert_midpoint_and_new_leaf(tree, diameter_path, midpoint_position, distance_to_midpoint, new_leaf_name, branch_length):
    parent = midpoint_position.up

    new_node = Tree(name="midpoint")
    new_node.dist = distance_to_midpoint
    midpoint_position.dist -= distance_to_midpoint

    if parent:
        parent.add_child(new_node)
    else:
        print("Error: Midpoint calculation attempts to modify the root.")
        return tree

    new_node.add_child(midpoint_position.detach())

    print(f"Inserting new leaf '{new_leaf_name}' with branch length {branch_length} at midpoint")
    new_leaf = Tree(name=new_leaf_name)
    new_leaf.dist = branch_length
    new_node.add_child(new_leaf)

    # Validate the distances from midpoint to both ends of the diameter path
    distance_to_leaf1 = tree.get_distance(new_node, diameter_path[0])
    distance_to_leaf2 = tree.get_distance(new_node, diameter_path[-1])
    print(f"Distance from midpoint to {diameter_path[0].name}: {distance_to_leaf1}")
    print(f"Distance from midpoint to {diameter_path[-1].name}: {distance_to_leaf2}")

    if abs(distance_to_leaf1 - half_distance) > 1e-6 or abs(distance_to_leaf2 - half_distance) > 1e-6:
        print(f"Warning: Midpoint distances are not equal. distance_to_leaf1: {distance_to_leaf1}, distance_to_leaf2: {distance_to_leaf2}, half_distance: {half_distance}")

    return tree

# Example
newick = "((A:1.5,B:1.2):1.5,(C:1.9,D:0.1):0.5);"
tree = Tree(newick, format=1)

# Label internal nodes
label_internal_nodes(tree)

print("Original tree with labeled internal nodes:")
print(tree.write(format=1))

temporary_leaves = {tree & "A", tree & "C", tree & "D"}

new_leaf_name = "New_leaf"
branch_length = 1.5

midpoint_position, diameter_path, half_distance, distance_to_midpoint = compute_midpoint(tree, temporary_leaves)
tree = insert_midpoint_and_new_leaf(tree, diameter_path, midpoint_position, distance_to_midpoint, new_leaf_name, branch_length)

print("Updated tree:")
print(tree.write(format=1))

