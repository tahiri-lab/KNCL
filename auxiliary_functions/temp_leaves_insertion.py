# Case 1: insert a temporary leaf on the terminal branch of the target leaf

from ete3 import Tree

def insert_leaf(newick, target_leaf, new_leaf, new_length, dist):
    # Parse the Newick format tree
    tree = Tree(newick, format=1)
    
    # Find the target leaf node
    node = tree.search_nodes(name=target_leaf)[0]
    
    # Ensure the distance is not greater than the branch length of the target leaf
    if dist > node.dist:
        raise ValueError("Distance is greater than the branch length of the target leaf.")
    
    # If dist is less than the node's distance, we need to create a new internal node
    if dist < node.dist:
        # Calculate the remaining distance of the original branch after insertion
        remaining_dist = node.dist - dist
        
        # Detach the node temporarily
        parent = node.up
        node.detach()
        
        # Create a new internal node and attach it to the original parent
        new_internal_node = parent.add_child(dist=dist)
        
        # Reattach the original node and add the new leaf
        node.dist = remaining_dist
        new_internal_node.add_child(node)
        new_internal_node.add_child(name=new_leaf, dist=new_length)
    else:
        # If dist equals the node's distance, attach directly to the node
        node.add_child(name=new_leaf, dist=new_length)

    # Return the modified tree in Newick format
    return tree

# Example
newick = "(A:0.5,(B:0.3,C:0.4):0.7);"
target_leaf = "C"
new_leaf = "L"
new_length = 0.2
dist = 0.1

# Insert the leaf and print the new tree
new_tree = insert_leaf(newick, target_leaf, new_leaf, new_length, dist)
print(new_tree)

# Case 2. The target distance is greater than the length of the terminal branch of the target leaf

from ete3 import Tree

def insert_leaf(newick, target_leaf, new_leaf, new_length, dist):
    tree = Tree(newick, format=1)
    
    # Find the target leaf node
    target_node = tree.search_nodes(name=target_leaf)[0]
    current_node = target_node
    accumulated_distance = 0

    # Trace back up the tree to find the correct insertion point
    while current_node.up:
        if accumulated_distance + current_node.dist > dist:
            # Calculate the exact point where the new node should be inserted
            insert_distance = dist - accumulated_distance
            break
        accumulated_distance += current_node.dist
        current_node = current_node.up

    if not current_node.up and accumulated_distance + current_node.dist < dist:
        raise ValueError("Specified distance exceeds the available path in the tree.")

    # Perform the insertion
    excess_length = current_node.dist - insert_distance
    new_internal_node = current_node.up.add_child(dist=insert_distance)
    
    current_node.detach()
    new_internal_node.add_child(current_node, dist=excess_length)
    new_internal_node.add_child(name=new_leaf, dist=new_length)

    return tree

# Example
newick = "(A:0.5,((B:0.3,C:1.9):0.7,D:0.5):0.8);"
target_leaf = "C"
new_leaf = "L"
new_length = 0.2
dist = 0.3  # Distance from the existing leaf C
dist2 = 2.0  # A larger distance to test more extensive tree traversal

# Insert the leaf and print the new tree
new_tree1 = insert_leaf(newick, target_leaf, new_leaf, new_length, dist)
new_tree2 = insert_leaf(newick, target_leaf, new_leaf, new_length, dist2)
print("Insertion at smaller distance:", new_tree1)
print(new_tree1.write(format=1))
print("Insertion at larger distance:", new_tree2)
print(new_tree2.write(format=1))
