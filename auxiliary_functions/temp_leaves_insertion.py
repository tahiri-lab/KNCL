# Temporary leaf insertion
# Make sure you have ete3 installed

from ete3 import Tree

def insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist):
    tree = Tree(newick, format=1)
    target_node = tree.search_nodes(name=target_leaf)[0]
    insertion_points = []  # Store nodes where inserts will be made
    visited_nodes = set()  # Set to track visited nodes

    def insert_leaf_at_node(parent_node, insert_distance, previous_node):
        if insert_distance == 0:
            new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
            parent_node.add_child(name=new_leaf_name, dist=round(new_length, 8))
            insertion_points.append(parent_node)
            visited_nodes.add(parent_node)
            visited_nodes.add(new_leaf_name)
        else:
            excess_length = round(parent_node.dist,8) - round(insert_distance, 8)
            new_internal_node = parent_node.up.add_child(dist=round(insert_distance, 8)) if parent_node.up else tree.add_child(dist=round(insert_distance, 8))
            parent_node.detach()
            new_internal_node.add_child(parent_node, dist=excess_length)
            new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=round(new_length, 8))
            insertion_points.append(new_internal_node)
            visited_nodes.add(new_internal_node)
            visited_nodes.add(new_internal_node.children[1])  # Newly added leaf node

    def dfs(node, accumulated_distance, previous_node=None, previous_distance=0):
        if node in visited_nodes:
            return
        visited_nodes.add(node)

        if round(accumulated_distance, 8) >= dist:
            insert_distance = round(accumulated_distance,8) - round(dist, 8)
            if insert_distance == 0:
                # Insert leaf directly at the current node
                insert_leaf_at_node(node, insert_distance, node)
            elif node.is_leaf():
                # Insert new leaf at the correct branch between the node and its parent
                insert_leaf_at_node(node, insert_distance, node)
            else:
                # Insert new leaf at the correct branch between the previous node and the current node
                insert_leaf_at_node(previous_node, round(previous_distance - insert_distance, 8), previous_node)
            return

        # Traverse children
        for child in node.children:
            if child not in visited_nodes:  # Avoid traversing newly added leaves
                dfs(child, accumulated_distance + child.dist, node, child.dist)

        # Traverse parent
        if node.up and node.up not in visited_nodes:  # Avoid traversing newly added internal nodes
            dfs(node.up, accumulated_distance + node.dist, node, node.dist)

    # Start the DFS traversal from the target node
    dfs(target_node, 0)

    # Output the final results
    if insertion_points:
        print("Final tree with all inserted leaves:")
        print(tree.write(format=1))
        print(tree)
    else:
        print("No valid insertion points were found based on the specified distance.")

# Example
newick = "(A:0.5,((B:0.3,C:1.9)Node1:0.7,D:0.5)Node2:0.8);"
target_leaf = "C"
new_leaf_base_name = "L"
new_length = 0.2
dist = 2.7  # Distance to test for multiple possible insertions

# Insert new leaves and check the tree structure
insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist)
