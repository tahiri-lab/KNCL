# Finding the midpoint among a set of temporary leaves TL and inserting a new leaf at the midpoint

from ete3 import Tree

def find_farthest_leaf(tree, start, temporary_leaves):
    max_distance = 0
    farthest_leaf = start
    for leaf in temporary_leaves:
        if leaf is not start:
            distance = tree.get_distance(start, leaf)
            if distance > max_distance:
                max_distance = distance
                farthest_leaf = leaf
    print(f"Farthest leaf from {start.name}: {farthest_leaf.name} at distance {max_distance}")
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

    # Find the common ancestor
    lca = None
    for n1 in path1:
        if n1 in path2:
            lca = n1
            break
    
    path = []
    for node in path1:
        path.append(node)
        if node == lca:
            break
    
    lca_index = path2.index(lca)
    path.extend(reversed(path2[:lca_index]))
    
    path_names = [n.name for n in path]
    print(f"Diameter path: {' -> '.join(path_names)}")
    return path

def label_internal_nodes(tree):
    counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.name = f"Node{counter}"
            counter += 1

def compute_midpoint(tree, temporary_leaves):
    start = next(iter(temporary_leaves))
    leaf1, dist1 = find_farthest_leaf(tree, start, temporary_leaves)
    leaf2, dist2 = find_farthest_leaf(tree, leaf1, temporary_leaves)
    path = find_path(leaf1, leaf2)
    total_distance = dist2
    half_distance = total_distance / 2

    cumulative_distance = 0
    for node in path:
        cumulative_distance += node.dist
        if cumulative_distance > half_distance:
            excess = cumulative_distance - half_distance
            print(f"Midpoint at {node.name}, excess: {excess}, node distance: {node.dist}")
            return node, excess, path

def insert_midpoint_and_new_leaf(tree, midpoint_node, excess, new_leaf_name, branch_length, path):
    print(f"Inserting new leaf at {midpoint_node.name} with excess {excess}")
    if excess == 0:
        new_leaf = Tree(name=new_leaf_name)
        new_leaf.dist = branch_length
        midpoint_node.add_child(new_leaf)
    else:
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
        #new_node.dist = excess
        new_node.dist = distance_to_midpoint

        parent = midpoint_node.up
        parent.remove_child(midpoint_node)
        parent.add_child(new_node)
        #midpoint_node.dist = midpoint_node.dist - excess
        midpoint_node.dist = distance_from_midpoint_to_leaf
        new_node.add_child(midpoint_node)

        new_leaf = Tree(name=new_leaf_name)
        new_leaf.dist = branch_length
        new_node.add_child(new_leaf)

    return tree

# Example
newick = "(A:2,(D:0.5,((C:1.9,(B:0.2,C1:0.2):0.1):0.1,C2:0.2):0.6):0.8);"
tree = Tree(newick, format=1)
label_internal_nodes(tree)
print("Original tree with labeled nodes:")
print(tree)

temporary_leaves = {tree & "C1", tree & "A"}
new_leaf_name = "New_leaf"
branch_length = 1.5

midpoint_node, excess, path = compute_midpoint(tree, temporary_leaves)
tree = insert_midpoint_and_new_leaf(tree, midpoint_node, excess, new_leaf_name, branch_length, path)

print("Updated tree:")
print(tree)
print(tree.write(format=1))

# Debugging part
from ete3 import Tree

def find_farthest_leaf(tree, start, temporary_leaves):
    max_distance = 0
    farthest_leaf = start
    for leaf in temporary_leaves:
        if leaf is not start:
            distance = tree.get_distance(start, leaf)
            if distance > max_distance:
                max_distance = distance
                farthest_leaf = leaf
    print(f"Farthest leaf from {start.name}: {farthest_leaf.name} at distance {max_distance}")
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

    # Find the common ancestor
    lca = None
    for n1 in path1:
        if n1 in path2:
            lca = n1
            break

    path = []
    for node in path1:
        path.append(node)
        if node == lca:
            break

    lca_index = path2.index(lca)
    path.extend(reversed(path2[:lca_index]))

    path_names = [n.name for n in path]
    print(f"Diameter path: {' -> '.join(path_names)}")
    return path

def label_internal_nodes(tree):
    counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            node.name = f"Node{counter}"
            counter += 1

def compute_midpoint(tree, temporary_leaves):
    start = next(iter(temporary_leaves))
    leaf1, dist1 = find_farthest_leaf(tree, start, temporary_leaves)
    leaf2, dist2 = find_farthest_leaf(tree, leaf1, temporary_leaves)

    # Determine the leaf closest to the root to use as the consistent endpoint
    root_distances = {leaf1: tree.get_distance(tree.get_tree_root(), leaf1), leaf2: tree.get_distance(tree.get_tree_root(), leaf2)}
    endpoint = max(root_distances, key=root_distances.get)

    if endpoint == leaf2:
        leaf1, leaf2 = leaf2, leaf1

    path = find_path(leaf1, leaf2)
    total_distance = dist2
    half_distance = total_distance / 2

    cumulative_distance = 0
    for i, node in enumerate(path):
        if i > 0:
            cumulative_distance += path[i - 1].dist
        print(f"Node: {node.name}, Dist: {node.dist}, Cumulative distance: {cumulative_distance}")
        if cumulative_distance >= half_distance:
            excess = cumulative_distance - half_distance
            prev_node = path[i - 1]
            print(f"Midpoint between {prev_node.name} and {node.name}, excess: {excess}")
            return prev_node, node, excess, path

def insert_midpoint_and_new_leaf(tree, prev_node, curr_node, excess, new_leaf_name, branch_length, path):
    print(f"Inserting new leaf between {prev_node.name} and {curr_node.name} with excess {excess}")

    # Calculate distances
    original_dist = prev_node.get_distance(curr_node)
    distance_to_midpoint = excess
    distance_from_midpoint_to_leaf = original_dist - excess

    print(f"Original distance: {original_dist}, Distance to midpoint: {distance_to_midpoint}, Distance from midpoint to leaf: {distance_from_midpoint_to_leaf}")

    # Insert new midpoint node
    new_node = Tree(name="midpoint")
    new_node.dist = distance_to_midpoint

    # Attach new midpoint node correctly
    if prev_node.up == curr_node:
        parent = curr_node
        child = prev_node
    else:
        parent = prev_node
        child = curr_node

    print(f"Parent node: {parent.name}, Child node: {child.name}")

    # Remove the child from the parent
    parent.remove_child(child)
    print(f"Removed child {child.name} from parent {parent.name}")

    # Add the new midpoint node to the parent
    parent.add_child(new_node)
    print(f"Added new midpoint node 'midpoint' to parent {parent.name}")

    # Reattach the child to the new midpoint node
    new_node.add_child(child)
    child.dist = distance_from_midpoint_to_leaf
    print(f"Reattached current node '{child.name}' to new midpoint node with distance {child.dist}")

    # Add the new leaf to the new midpoint node
    new_leaf = Tree(name=new_leaf_name)
    new_leaf.dist = branch_length
    new_node.add_child(new_leaf)
    print(f"Added new leaf '{new_leaf_name}' to midpoint node 'midpoint' with distance {new_leaf.dist}")

    print(f"Inserted new internal node 'midpoint' with new leaf '{new_leaf.name}'")

    return tree

# Example
newick = "((A:2,A2:1.3):0.85,(D:0.5,((C:1.9,(B:0.2,C1:0.2):0.1):0.1,C2:0.2):0.6):0.8);"
tree = Tree(newick, format=1)
label_internal_nodes(tree)
print("Original tree with labeled nodes:")
print(tree)

temporary_leaves = {tree & "C2", tree & "A"}
new_leaf_name = "New_leaf"
branch_length = 1.5

prev_node, curr_node, excess, path = compute_midpoint(tree, temporary_leaves)
tree = insert_midpoint_and_new_leaf(tree, prev_node, curr_node, excess, new_leaf_name, branch_length, path)

print("Updated tree:")
print(tree)
print(tree.write(format=1))
