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
    new_internal_node = current_node.up.add_child(dist=excess_length)
    
    current_node.detach()
    new_internal_node.add_child(current_node, dist=insert_distance)
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

Case 3. Universal approach

from ete3 import Tree

def insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist):
    tree = Tree(newick, format=1)
    target_node = tree.search_nodes(name=target_leaf)[0]
    insertion_points = []  # Store nodes and distances where inserts will be made
    ignored_nodes = set()  # Track nodes to be ignored during traversal

    # Function to perform the insertion at a given node with calculated distances
    def perform_insertion(node, insert_distance):
        excess_length = node.dist - insert_distance
        if node.up:
            new_internal_node = node.up.add_child(dist=insert_distance)
            node.detach()
            new_internal_node.add_child(node, dist=excess_length)
            new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
            new_leaf = new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_internal_node)  # Record the insertion
            ignored_nodes.add(new_internal_node)  # Ignore this newly created node in further traversals
            ignored_nodes.add(new_leaf)  # Also ignore the new leaf
            print(f"Inserted '{new_leaf_name}' at node '{node.name}' with insert distance {insert_distance} and excess length {excess_length}")
            return new_internal_node
        return None

    # Traverse upwards from the target node
    def traverse_upwards(node, accumulated_distance):
        current_distance = accumulated_distance + node.dist
        print(f"Traversing upwards through '{node.name}' with current distance {current_distance}")
        if current_distance >= dist:
            insert_distance = current_distance - dist
            if insert_distance <= node.dist:
                inserted_node = perform_insertion(node, insert_distance)
                if inserted_node:
                    return True  # Stop traversing upwards once a valid insertion is made
        if node.up:
            return traverse_upwards(node.up, current_distance)
        return False

    # Traverse downwards from the target node
    def traverse_downwards(node, accumulated_distance):
        if node in ignored_nodes:
            return  # Skip nodes that should be ignored

        current_distance = accumulated_distance + node.dist
        print(f"Traversing downwards from '{node.name}' with current distance {current_distance}")
        if current_distance >= dist:
            insert_distance = current_distance - dist
            if insert_distance <= node.dist:
                inserted_node = perform_insertion(node, insert_distance)
                if inserted_node:
                    ignored_nodes.add(inserted_node)
        for child in node.children:
            if child not in ignored_nodes:
                traverse_downwards(child, current_distance)  # Accumulate the distance properly

    # Start traversal from the target node upwards
    traverse_upwards(target_node, 0)

    # Start traversal from the parent node of the target node downwards
    if target_node.up:
        for sibling in target_node.up.children:
            if sibling != target_node:
                traverse_downwards(sibling, target_node.dist)

    # Output the results
    if insertion_points:
        print("Final tree with all inserted leaves:")
        print(tree.write(format=1))
        print(tree)
    else:
        print("No valid insertion points were found based on the specified distance.")

# Example
newick = "(A:0.5,((B:0.3,C:1.9):0.7,D:0.5):0.8);"
target_leaf = "C"
new_leaf_base_name = "L"
new_length = 0.2
dist = 2.1  # Distance to test for multiple possible insertions

# Insert new leaves and check the tree structure
insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist)
