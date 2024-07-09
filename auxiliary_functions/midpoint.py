# Finding the midpoint among a set of temporary leaves TL

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
    return midpoint_position, diameter_path

def insert_midpoint_and_new_leaf(tree, diameter_path, new_leaf_name, branch_length):
    total_distance = sum(node.dist for node in diameter_path)
    half_distance = total_distance / 2
    cumulative_distance = 0
    for i, node in enumerate(diameter_path):
        cumulative_distance += node.dist
        if cumulative_distance >= half_distance:
            break

    midpoint_position = diameter_path[i]
    parent = midpoint_position.up

    new_node = Tree(name="midpoint")
    new_node.dist = midpoint_position.dist / 2
    midpoint_position.dist /= 2

    parent.add_child(new_node)
    new_node.add_child(midpoint_position.detach())

    print(f"Inserting new leaf '{new_leaf_name}' with branch length {branch_length} at midpoint")
    new_leaf = Tree(name=new_leaf_name)
    new_leaf.dist = branch_length
    new_node.add_child(new_leaf)

    return tree

# Example
if __name__ == "__main__":
    newick = "((A:1.5,B:1.2):2,(C:1.9,D:0.1):2);"
    tree = Tree(newick, format=1)

    print("Original tree:")
    print(tree)

    temporary_leaves = {tree & "A", tree & "C", tree & "D"}
    new_leaf_name = "New_leaf"
    branch_length = 1.5

    midpoint_position, diameter_path = compute_midpoint(tree, temporary_leaves)
    tree = insert_midpoint_and_new_leaf(tree, diameter_path, new_leaf_name, branch_length)

    print("Updated tree:")
    print(tree)
    print(tree.write(format=1))
