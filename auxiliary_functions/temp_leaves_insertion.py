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
dist = 1.9  # Distance to test for multiple possible insertions

# Insert new leaves and check the tree structure
insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist)

# Debugging

from ete3 import Tree

def insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist, tolerance=1e-10):
    tree = Tree(newick, format=1)
    target_node = tree.search_nodes(name=target_leaf)[0]
    insertion_points = []  # Store nodes where inserts will be made
    visited_nodes = set()  # Set to track visited nodes

    # Label internal nodes for easier tracking
    internal_node_counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and not node.name:
            node.name = f"Node{internal_node_counter}"
            internal_node_counter += 1

    def insert_leaf_at_node(parent_node, insert_distance, previous_node):
        print(f"Attempting to insert at node '{parent_node.name}' with insert distance {insert_distance}, coming from '{previous_node.name if previous_node else 'None'}'")
        if abs(insert_distance) < tolerance:
            new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
            parent_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(parent_node)
            visited_nodes.add(parent_node)
            visited_nodes.add(new_leaf_name)
            print(f"Inserted '{new_leaf_name}' at existing node '{parent_node.name}' with no new internal node.")
        else:
            excess_length = parent_node.dist - insert_distance
            if excess_length < 0:
                print(f"Invalid insertion leading to negative branch length: {excess_length}")
                return False

            if parent_node.up:
                new_internal_node = parent_node.up.add_child(dist=insert_distance)
                parent_node.detach()
                new_internal_node.add_child(parent_node, dist=excess_length)
                new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
                new_internal_node.add_child(name=new_leaf_name, dist=new_length)
                insertion_points.append(new_internal_node)
                visited_nodes.add(new_internal_node)
                visited_nodes.add(new_internal_node.children[1])  # Newly added leaf node
                print(f"Inserted '{new_leaf_name}' between '{parent_node.name}' and '{parent_node.up.name}' with insert distance {insert_distance} and excess length {excess_length}")
            else:
                # Handling the case when the parent node is the root
                new_internal_node = tree.add_child(dist=insert_distance)
                parent_node.detach()
                new_internal_node.add_child(parent_node, dist=excess_length)
                new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
                new_internal_node.add_child(name=new_leaf_name, dist=new_length)
                insertion_points.append(new_internal_node)
                visited_nodes.add(new_internal_node)
                visited_nodes.add(new_internal_node.children[1])  # Newly added leaf node
                print(f"Inserted '{new_leaf_name}' at root between '{parent_node.name}' and the root with insert distance {insert_distance} and excess length {excess_length}")
            return True

    def bfs(node, accumulated_distance):
        queue = [(node, accumulated_distance, None, 0, [])]
        while queue:
            current_node, current_dist, prev_node, prev_dist, path = queue.pop(0)
            if current_node in visited_nodes:
                continue
            visited_nodes.add(current_node)
            current_path = path + [current_node.name]

            print(f"Traversing '{current_node.name}' with accumulated distance: {current_dist}. Path: {' -> '.join(current_path)}")
            if round(current_dist, 8) >= dist:
                insert_distance = round(current_dist, 8) - round(dist, 8)
                if abs(insert_distance) < tolerance:
                    insert_distance = 0
                if insert_distance == 0:
                    # Insert leaf directly at the current node
                    if not insert_leaf_at_node(current_node, insert_distance, prev_node):
                        return
                elif current_node.is_leaf():
                    # Insert new leaf at the correct branch between the node and its parent
                    if not insert_leaf_at_node(current_node, insert_distance, prev_node):
                        return
                else:
                    # Insert new leaf at the correct branch between the previous node and the current node
                    print(f"Checking insertion between previous node '{prev_node.name if prev_node else 'None'}' and current node '{current_node.name}' with distances {prev_dist} - {insert_distance}")
                    if not insert_leaf_at_node(prev_node, prev_dist - insert_distance, current_node):
                        return
                continue

            # Traverse children
            for child in current_node.children:
                if child not in visited_nodes:  # Avoid traversing newly added leaves
                    queue.append((child, current_dist + child.dist, current_node, child.dist, current_path))

            # Traverse parent
            if current_node.up and current_node.up not in visited_nodes:  # Avoid traversing newly added internal nodes
                queue.append((current_node.up, current_dist + current_node.dist, current_node, current_node.dist, current_path))

    def bfs_reverse(node, accumulated_distance):
        queue = [(node, accumulated_distance, None, 0, [])]
        while queue:
            current_node, current_dist, prev_node, prev_dist, path = queue.pop(0)
            if current_node in visited_nodes:
                continue
            visited_nodes.add(current_node)
            current_path = path + [current_node.name]

            print(f"Traversing '{current_node.name}' with accumulated distance: {current_dist}. Path: {' -> '.join(current_path)}")
            if round(current_dist, 8) >= dist:
                insert_distance = round(current_dist, 8) - round(dist, 8)
                if abs(insert_distance) < tolerance:
                    insert_distance = 0
                if insert_distance == 0:
                    # Insert leaf directly at the current node
                    if not insert_leaf_at_node(current_node, insert_distance, prev_node):
                        return
                elif current_node.is_leaf():
                    # Insert new leaf at the correct branch between the node and its parent
                    if not insert_leaf_at_node(current_node, insert_distance, prev_node):
                        return
                else:
                    # Insert new leaf at the correct branch between the previous node and the current node
                    print(f"Checking insertion between previous node '{prev_node.name if prev_node else 'None'}' and current node '{current_node.name}' with distances {prev_dist} - {insert_distance}")
                    if not insert_leaf_at_node(current_node, insert_distance, prev_node):
                        return
                continue

            # Traverse parent
            if current_node.up and current_node.up not in visited_nodes:  # Avoid traversing newly added internal nodes
                queue.append((current_node.up, current_dist + current_node.dist, current_node, current_node.dist, current_path))

            # Traverse children
            for child in current_node.children:
                if child not in visited_nodes:  # Avoid traversing newly added leaves
                    queue.append((child, current_dist + child.dist, current_node, child.dist, current_path))

    # Direct insertion case when dist is less than the terminal branch length of the target leaf
    if dist <= target_node.dist:
        print(f"Direct insertion at target leaf '{target_leaf}' with distance {dist}")
        insert_leaf_at_node(target_node, dist, target_node)

    # Use BFS for inward traversal and BFS reverse for outward traversal
    elif tree.get_common_ancestor(tree.get_leaves()) == target_node:
        bfs(target_node, 0)
    else:
        bfs_reverse(target_node, 0)

    # Round the final output distances for better readability
    def round_tree_distances(tree_node, decimals=8):
        for node in tree_node.traverse():
            node.dist = round(node.dist, decimals)

    round_tree_distances(tree)

    # Output the final results
    if insertion_points:
        print("Final tree with all inserted leaves:")
        print(tree.write(format=1))
        print(tree)
    else:
        print("No valid insertion points were found based on the specified distance.")

# First example
newick = "(((A:1.587,(F:1.110,(M:1.343,R:1.369):0.846):0.487):1.981,D:0.356):2.121,(B:1.936,(C:0.915,Q:1.201):2.101):0.912);"
target_leaf = "D"
new_leaf_base_name = "E"
new_length = 0.279
dist = 2.695936081694403

print(Tree(newick, format=1))
# Insert new leaves and check the tree structure
insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist)
print(Tree(newick, format=1))

# Second example
newick = "(((A:1.587,(F:1.110,(M:1.343,R:1.369):0.846):0.487):1.981,D:0.356):2.121,(B:1.936,(C:0.915,Q:1.201):2.101):0.912);"
target_leaf = "Q"
new_leaf_base_name = "E"
new_length = 0.279
dist = 3.0597060866386405

print(Tree(newick, format=1))
# Insert new leaves and check the tree structure
insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist)
print(Tree(newick, format=1))

# Third example
newick = "(((A:1.587,(F:1.110,(M:1.343,R:1.369):0.846):0.487):1.981,D:0.356):2.121,(B:1.936,(C:0.915,Q:1.201):2.101):0.912);"
target_leaf = "C"
new_leaf_base_name = "E"
new_length = 0.279
dist = 0.67138071756504

print(Tree(newick, format=1))
# Insert new leaves and check the tree structure
insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist)
print(Tree(newick, format=1))
