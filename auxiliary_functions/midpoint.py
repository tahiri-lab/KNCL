# Finding the midpoint among a set of temporary leaves TL

from ete3 import Tree

def bfs(start, temporary_leaves):
    from collections import deque

    queue = deque([start])
    distances = {start: 0}
    farthest_leaf = start

    while queue:
        current = queue.popleft()
        for child in current.get_children():
            if child in temporary_leaves and child not in distances:
                distances[child] = distances[current] + 1
                queue.append(child)
                if distances[child] > distances[farthest_leaf]:
                    farthest_leaf = child

    return farthest_leaf, distances

def find_path(distances, leaf1, leaf2):
    path = []
    current = leaf2
    while current != leaf1:
        path.append(current)
        current = current.up
    path.append(leaf1)
    path.reverse()
    return path

def compute_midpoint(tree, temporary_leaves):
    start = next(iter(temporary_leaves))
    leaf1, distances = bfs(start, temporary_leaves)
    leaf2, distances = bfs(leaf1, temporary_leaves)
    diameter_path = find_path(distances, leaf1, leaf2)
    half_distance = len(diameter_path) // 2
    midpoint_position = diameter_path[half_distance]

    return midpoint_position, diameter_path

def insert_midpoint(tree, diameter_path):
    half_distance = len(diameter_path) // 2
    midpoint_position = diameter_path[half_distance]
    parent = midpoint_position.up

    new_node = Tree(name="midpoint")
    new_node.dist = midpoint_position.dist / 2
    midpoint_position.dist /= 2

    parent.add_child(new_node)
    new_node.add_child(midpoint_position.detach())

    return tree

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

newick = "((A:1.5,B:1.2):2,(C:1.9,D:0.1):2);"
tree = Tree(newick, format=1)

temporary_leaves = {tree & "A", tree & "C", tree & "D"}

midpoint_position, diameter_path = compute_midpoint(tree, temporary_leaves)
tree = insert_midpoint(tree, diameter_path)

print(tree.write(format=1))
