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

def insert_leaf_extended(newick, target_leaf, new_leaf, new_length, dist):
    tree = Tree(newick, format=1)
    
    # Find the target leaf node
    target_node = tree.search_nodes(name=target_leaf)[0]

    # Traverse the tree to find a suitable insertion point
    # This function will attempt to move upward and then to distant leaves if necessary
    def traverse_for_insertion(node, current_dist, path_visited):
        if current_dist >= dist:
            # Calculate the exact insertion point on this branch
            parent_dist = node.up.dist if node.up else 0
            insert_dist = dist - (current_dist - parent_dist)
            new_internal_node = node.up.add_child(name="internal_node", dist=insert_dist)
            node.up.dist = parent_dist - insert_dist
            new_internal_node.add_child(name=new_leaf, dist=new_length)
            return True
        else:
            # Mark this node as visited
            path_visited.add(node)
            # Check all connections (parent and siblings)
            connections = [node.up] if node.up and node.up not in path_visited else []
            connections.extend([sib for sib in node.get_sisters() if sib not in path_visited])
            for conn in connections:
                if traverse_for_insertion(conn, current_dist + conn.dist, path_visited):
                    return True
            # Mark this node as unvisited before going back
            path_visited.remove(node)
            return False

    # Start the search from the target node, considering upward and sibling pathways
    if not traverse_for_insertion(target_node, target_node.dist, set()):
        raise ValueError("Specified distance exceeds available pathways from the target leaf.")

    return tree

# Example
newick = "(A:0.5,(B:0.3,C:0.4):0.7,(E:1.2,F:1.3):0.8);"
target_leaf = "C"
new_leaf = "L"
new_length = 0.2
dist = 0.3  # Distance from the existing leaf C

# Insert the leaf and print the new tree
new_tree = insert_leaf_extended(newick, target_leaf, new_leaf, new_length, dist)
print(new_tree)
