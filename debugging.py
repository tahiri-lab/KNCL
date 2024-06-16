#Debugging of the updated version of the k-NCL algorithm

from ete3 import Tree #do not forget to install ete3 if needed
import copy

def find_sd(tree, distinct_leaves):
    visited_nodes = set()
    subtree_roots = set()
    sd = []
    
    leaf_to_node = {leaf.name: leaf for leaf in tree.iter_leaves()}
    node_to_parent = {}

    def map_parents(node, parent=None):
        node_to_parent[node] = parent
        for child in node.children:
            map_parents(child, node)

    map_parents(tree)
    
    for leaf in distinct_leaves:
        node = leaf_to_node[leaf]
        current_node = node
        subtree_root = None
        
        while current_node and current_node not in visited_nodes:
            visited_nodes.add(current_node)
            leaves_current_node = {leaf.name for leaf in current_node.iter_leaves()}
            if leaves_current_node.issubset(distinct_leaves):
                subtree_root = current_node
            current_node = node_to_parent.get(current_node, None)
            
        if subtree_root and subtree_root not in subtree_roots:
            subtree_roots.add(subtree_root)
            sd.append(subtree_root)
    
    return sd

def compute_distance(tree, leaf1, leaf2):
    nodes1 = tree.search_nodes(name=leaf1)
    nodes2 = tree.search_nodes(name=leaf2)
    if not nodes1 or not nodes2:
        raise ValueError(f"One of the nodes {leaf1} or {leaf2} not found in the tree.")
    node1 = nodes1[0]
    node2 = nodes2[0]
    return node1.get_distance(node2)

def compute_adjustment_rate(tree1, tree2, common_leaves):
    common_leaves = list(common_leaves)
    sum_distances_tree1 = sum(compute_distance(tree1, l1, l2) for i, l1 in enumerate(common_leaves[:-1]) for l2 in common_leaves[i+1:])
    sum_distances_tree2 = sum(compute_distance(tree2, l1, l2) for i, l1 in enumerate(common_leaves[:-1]) for l2 in common_leaves[i+1:])
    
    if sum_distances_tree2 == 0:
        return 1.0  # Default adjustment rate if sum_distances_tree2 is zero, but usually it's not the case
    
    return sum_distances_tree1 / sum_distances_tree2

def k_nearest_common_leaves(tree, element, common_leaves, k):
    distances = [(leaf, compute_distance(tree, element, leaf)) for leaf in common_leaves]
    distances.sort(key=lambda x: x[1])
    nearest_leaves = [leaf for leaf, distance in distances[:k]]
    return nearest_leaves

def compute_midpoint(tree, temp_leaves):
    if len(temp_leaves) < 2:
        raise ValueError("Not enough temporary leaves to compute midpoint.")
    
    max_distance = 0
    distant_pair = (None, None)
    
    for i in range(len(temp_leaves) - 1):
        for j in range(i + 1, len(temp_leaves)):
            distance = compute_distance(tree, temp_leaves[i], temp_leaves[j])
            if distance > max_distance:
                max_distance = distance
                distant_pair = (temp_leaves[i], temp_leaves[j])
    
    if not distant_pair[0] or not distant_pair[1]:
        raise ValueError(f"Could not find distant pair. temp_leaves: {temp_leaves}, distant_pair: {distant_pair}")
    
    half_distance = max_distance / 2
    midpoint = traverse_tree(tree, distant_pair, half_distance)
    return midpoint

def traverse_tree(tree, distant_pair, half_distance):
    nodes1 = tree.search_nodes(name=distant_pair[0])
    nodes2 = tree.search_nodes(name=distant_pair[1])
    if not nodes1 or not nodes2:
        raise ValueError(f"One of the nodes {distant_pair[0]} or {distant_pair[1]} not found in the tree.")
    node1 = nodes1[0]
    node2 = nodes2[0]
    common_ancestor = node1.get_common_ancestor(node2)
    
    accumulated_distance = 0
    for node in common_ancestor.iter_descendants("preorder"):
        if node.dist:
            accumulated_distance += node.dist
            if accumulated_distance >= half_distance:
                return node

    raise ValueError("Could not find the midpoint in the tree.")

def insert_temp_leaves(tree, common_leaf, temp_distance):
    temp_leaves = []
    nodes = tree.search_nodes(name=common_leaf)
    if not nodes:
        raise ValueError(f"Node {common_leaf} not found in the tree.")
    common_leaf_node = nodes[0]
    path = common_leaf_node.get_ancestors()
    accumulated_distance = 0
    
    print(f"Inserting temporary leaf for {common_leaf} with temp_distance: {temp_distance}")
    
    for node in path:
        print(f"Traversing node: {node.name}, accumulated_distance: {accumulated_distance}, node.dist: {node.dist}")
        if node.dist:
            accumulated_distance += node.dist
            if accumulated_distance >= temp_distance:
                temp_leaf = Tree(name=f"temp_{common_leaf}")
                temp_leaf.dist = temp_distance
                temp_leaves.append(temp_leaf)
                node.add_child(temp_leaf)
                break
    
    if not temp_leaves:
        total_distance = sum(node.dist for node in path if node.dist)
        print(f"Total available distance: {total_distance}")
        if total_distance >= temp_distance:
            raise ValueError(f"Could not insert temporary leaf for {common_leaf} at distance {temp_distance}.")
        else:
            # Fallback: insert at the last possible position
            temp_leaf = Tree(name=f"temp_{common_leaf}")
            temp_leaf.dist = total_distance
            temp_leaves.append(temp_leaf)
            common_leaf_node.add_child(temp_leaf)
            print(f"Fallback insertion at total_distance: {total_distance}")
    
    return temp_leaves

def k_ncl_algorithm(tree1, tree2, k):
    common_leaves = set(leaf.name for leaf in tree1.iter_leaves()) & set(leaf.name for leaf in tree2.iter_leaves())
    completed_tree1 = tree1.copy()
    completed_tree2 = tree2.copy()
    
    for tree, other_tree in [(tree1, tree2), (tree2, tree1)]:
        distinct_leaves = set(leaf.name for leaf in tree.iter_leaves()) - common_leaves
        sd_elements = find_sd(tree, distinct_leaves)
        adjustment_rate = compute_adjustment_rate(tree, other_tree, common_leaves)
        
        for element in sd_elements:
            nearest_common_leaves = k_nearest_common_leaves(tree, element.name, common_leaves, k)
            temp_leaves = []
            
            for common_leaf in nearest_common_leaves:
                leaf_adjustment_rate = compute_adjustment_rate(tree, other_tree, [common_leaf])
                temp_distance = (compute_distance(tree, common_leaf, element.name) - element.dist) * leaf_adjustment_rate
                temp_leaves.extend(insert_temp_leaves(other_tree, common_leaf, temp_distance))
            
            print(f"Temp leaves for element {element.name}: {temp_leaves}")
            midpoint = compute_midpoint(other_tree, [leaf.name for leaf in temp_leaves])
            new_element = Tree(name=element.name)
            new_element.dist = element.dist * adjustment_rate
            midpoint.add_child(new_element)
    
    return completed_tree1, completed_tree2

# Testing with trees

tree1 = Tree("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6,E:0.7);")
tree2 = Tree("((A:0.15,(F:0.25,G:0.35):0.45):0.55,(C:0.65,D:0.75):0.85);")

completed_tree1, completed_tree2 = k_ncl_algorithm(tree1, tree2, k=3)
completed_tree1.write(outfile="completed_tree1.nwk")
completed_tree2.write(outfile="completed_tree2.nwk")

# Print the completed trees to check the output
print(completed_tree1)
print(completed_tree2)
