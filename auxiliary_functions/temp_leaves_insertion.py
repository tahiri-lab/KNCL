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
