# Temporary leaf insertion
# Make sure you have ete3 installed

from ete3 import Tree

def insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist, tolerance=1e-10):
    tree = Tree(newick, format=1)
    target_node = tree.search_nodes(name=target_leaf)[0]
    insertion_points = []
    visited_nodes = set()

    # Label internal nodes and branches for easier tracking
    internal_node_counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and not node.name:
            node.name = f"Node{internal_node_counter}"
            internal_node_counter += 1

    def robust_insert_leaf_at_node(current_node, insert_distance, previous_node, original_branch_distance, toward_root=False):
        excess_length = original_branch_distance - insert_distance

        if excess_length < 0:
            excess_length = 0

        # Handle traversal toward the root by ensuring correct branch selection
        if toward_root:
            temp = current_node
            current_node = previous_node
            previous_node = temp

        # Detach the previous node from its parent (the node leading to the root)
        parent = previous_node.up
        if parent:
            previous_node.detach()

        if parent is None:
            new_internal_node = tree.add_child(dist=excess_length)
            current_node.detach()
            new_internal_node.add_child(current_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_internal_node)
            visited_nodes.add(new_internal_node)
            visited_nodes.add(new_leaf_name)
        else:
            new_internal_node = parent.add_child(dist=excess_length)
            new_internal_node.add_child(previous_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_internal_node)
            visited_nodes.add(new_internal_node)

        # Post-insertion validation
        correct_insertion = validate_insertion_path(current_node, new_internal_node, previous_node, original_branch_distance)
        if not correct_insertion:
            print(f"Error: Insertion point verification failed between '{previous_node.name}' and '{current_node.name}'")
        return correct_insertion

    def insert_leaf_at_terminal(current_node, insert_distance):
        excess_length = current_node.dist - insert_distance
        if excess_length < 0:
            excess_length = 0

        parent = current_node.up
        if parent:
            current_node.detach()
            new_internal_node = parent.add_child(dist=excess_length)
            new_internal_node.add_child(current_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_internal_node)
            visited_nodes.add(new_internal_node)
            visited_nodes.add(new_leaf_name)
        else:
            return False

        return True

    def bfs(node, accumulated_distance):
        queue = [(node, accumulated_distance, None, 0, [], False)]
        while queue:
            current_node, current_dist, prev_node, prev_dist, path, toward_root = queue.pop(0)
            if current_node in visited_nodes:
                continue
            visited_nodes.add(current_node)
            current_path = path + [current_node.name]

            if round(current_dist, 8) >= dist:
                insert_distance = round(current_dist, 8) - round(dist, 8)
                if abs(insert_distance) < tolerance:
                    insert_distance = 0
                if insert_distance == 0:
                    if not robust_insert_leaf_at_node(current_node, insert_distance, prev_node, current_node.dist, toward_root):
                        return
                elif current_node.is_leaf():
                    if not insert_leaf_at_terminal(current_node, insert_distance):
                        return
                else:
                    if not robust_insert_leaf_at_node(prev_node, prev_dist - insert_distance, current_node, prev_dist, toward_root):
                        return
                continue

            for child in current_node.children:
                if child not in visited_nodes:
                    queue.append((child, current_dist + child.dist, current_node, child.dist, current_path, False))

            if current_node.up and current_node.up not in visited_nodes:
                queue.append((current_node.up, current_dist + current_node.dist, current_node, current_node.dist, current_path, True))

    def validate_insertion_path(current_node, new_internal_node, previous_node, original_branch_distance):
        # Verifies if the insertion happened between the correct nodes
        distance_check = current_node.get_distance(new_internal_node) + new_internal_node.get_distance(previous_node)
        return abs(distance_check - original_branch_distance) < tolerance

    if dist <= target_node.dist:
        insert_leaf_at_terminal(target_node, dist)
    else:
        bfs(target_node, 0)

    if insertion_points:
        print(tree.write(format=1))
        print(tree)
    else:
        print("No valid insertion points were found based on the specified distance.")

# Example
newick = "(((A:1.587,(F:1.110,(M:1.343,R:1.369):0.846):0.487):1.981,D:0.356):2.121,(B:1.936,(C:0.915,Q:1.201):2.101):0.912);"
target_leaf = "Q"
new_leaf_base_name = "temp"
new_length = 0.279
dist = 3.0597060866386405

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
                new_internal_node = parent_node.up.add_child(dist=excess_length)
                parent_node.detach()
                new_internal_node.add_child(parent_node, dist=insert_distance)
                new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
                new_internal_node.add_child(name=new_leaf_name, dist=new_length)
                insertion_points.append(new_internal_node)
                visited_nodes.add(new_internal_node)
                visited_nodes.add(new_internal_node.children[1])  # Newly added leaf node
                # Swap the assigned lengths
                print(f"Inserted '{new_leaf_name}' between '{parent_node.name}' and '{parent_node.up.name}' with insert distance {insert_distance} and excess length {excess_length}")
            else:
                # Handling the case when the parent node is the root
                new_internal_node = tree.add_child(dist=excess_length)
                parent_node.detach()
                new_internal_node.add_child(parent_node, dist=insert_distance)
                new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
                new_internal_node.add_child(name=new_leaf_name, dist=new_length)
                insertion_points.append(new_internal_node)
                visited_nodes.add(new_internal_node)
                visited_nodes.add(new_internal_node.children[1])  # Newly added leaf node
                print(f"Inserted '{new_leaf_name}' at root between '{parent_node.name}' and the root with insert distance {insert_distance} and excess length {excess_length}")
            return True

    def bfs_case_1(node, accumulated_distance):
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

            # Traverse children
            for child in current_node.children:
                if child not in visited_nodes:  # Avoid traversing newly added leaves
                    queue.append((child, current_dist + child.dist, current_node, child.dist, current_path))

            # Traverse parent
            if current_node.up and current_node.up not in visited_nodes:  # Avoid traversing newly added internal nodes
                queue.append((current_node.up, current_dist + current_node.dist, current_node, current_node.dist, current_path))

    def bfs_case_2(node, accumulated_distance):
        # BFS traversal with a focus on trees where the leaf is on a different side of the root
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

    # Direct insertion case when dist is less than the terminal branch length of the target leaf
    if dist <= target_node.dist:
        print(f"Direct insertion at target leaf '{target_leaf}' with distance {dist}")
        insert_leaf_at_node(target_node, dist, target_node)
    else:
        # Determine which BFS method to use based on tree topology
        if tree.get_distance(tree.get_tree_root(), target_node) < tree.get_farthest_leaf()[1] / 2:
            bfs_case_1(target_node, 0)
        else:
            bfs_case_2(target_node, 0)

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

# Example
newick = "(((A:1.587,(F:1.110,(M:1.343,R:1.369):0.846):0.487):1.981,D:0.356):2.121,(B:1.936,(C:0.915,Q:1.201):2.101):0.912);"
target_leaf = "D"
new_leaf_base_name = "E"
new_length = 0.279
dist = 2.695936081694403

print(Tree(newick, format=1))
# Insert new leaves and check the tree structure
insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist)
print(Tree(newick, format=1))

# Updated traversal case
from ete3 import Tree

def insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist, tolerance=1e-10):
    tree = Tree(newick, format=1)
    target_node = tree.search_nodes(name=target_leaf)[0]
    insertion_points = []
    visited_nodes = set()

    # Label internal nodes and branches for easier tracking
    internal_node_counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and not node.name:
            node.name = f"Node{internal_node_counter}"
            internal_node_counter += 1

    def robust_insert_leaf_at_node(current_node, insert_distance, previous_node, original_branch_distance, toward_root=False):
        print(f"\nAttempting insertion between nodes:")
        print(f"Current node: {current_node.name}")
        print(f"Previous node: {previous_node.name}")
        print(f"Original branch distance: {original_branch_distance}")
        print(f"Insertion distance: {insert_distance}")

        excess_length = original_branch_distance - insert_distance
        print(f"Calculated excess length: {excess_length}")

        if excess_length < 0:
            excess_length = 0

        # Handle traversal toward the root by ensuring correct branch selection
        if toward_root:
            print("Handling traversal toward the root...")
            temp = current_node
            current_node = previous_node
            previous_node = temp

        # Detach the previous node from its parent (the node leading to the root)
        parent = previous_node.up
        if parent:
            previous_node.detach()

        if parent is None:
            print("Handling root case")
            new_internal_node = tree.add_child(dist=excess_length)
            current_node.detach()
            new_internal_node.add_child(current_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_internal_node)
            visited_nodes.add(new_internal_node)
            visited_nodes.add(new_leaf_name)
        else:
            print(f"Normal case: Adding new internal node between '{previous_node.name}' and its parent.")
            #new_internal_node = parent.add_child(dist=insert_distance)
            #new_internal_node.add_child(previous_node, dist=excess_length)
            new_internal_node = parent.add_child(dist=excess_length)
            new_internal_node.add_child(previous_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_internal_node)
            visited_nodes.add(new_internal_node)
            print(f"Inserted leaf '{new_leaf_name}' between '{previous_node.name}' and '{current_node.name}'")

        # Post-insertion validation
        correct_insertion = validate_insertion_path(current_node, new_internal_node, previous_node, original_branch_distance)
        if not correct_insertion:
            print(f"Error: Insertion point verification failed between '{previous_node.name}' and '{current_node.name}'")
        return correct_insertion

    def insert_leaf_at_terminal(current_node, insert_distance):
        print(f"\nInserting at terminal node '{current_node.name}' with insert distance {insert_distance}")
        excess_length = current_node.dist - insert_distance
        print(f"Terminal node excess length: {excess_length}")
        if excess_length < 0:
            excess_length = 0

        parent = current_node.up
        if parent:
            current_node.detach()
            new_internal_node = parent.add_child(dist=excess_length)
            new_internal_node.add_child(current_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_internal_node)
            visited_nodes.add(new_internal_node)
            visited_nodes.add(new_leaf_name)
            print(f"Inserted '{new_leaf_name}' at terminal node '{current_node.name}' with insert distance {insert_distance} and excess length {excess_length}")
        else:
            print("Unexpected case: trying to insert at terminal root leaf.")
            return False

        return True

    def bfs(node, accumulated_distance):
        queue = [(node, accumulated_distance, None, 0, [], False)]
        while queue:
            current_node, current_dist, prev_node, prev_dist, path, toward_root = queue.pop(0)
            if current_node in visited_nodes:
                continue
            visited_nodes.add(current_node)
            current_path = path + [current_node.name]

            print(f"\nTraversing '{current_node.name}' with accumulated distance: {current_dist}. Path: {' -> '.join(current_path)}")
            if round(current_dist, 8) >= dist:
                insert_distance = round(current_dist, 8) - round(dist, 8)
                if abs(insert_distance) < tolerance:
                    insert_distance = 0
                if insert_distance == 0:
                    print("Direct insertion scenario triggered")
                    if not robust_insert_leaf_at_node(current_node, insert_distance, prev_node, current_node.dist, toward_root):
                        return
                elif current_node.is_leaf():
                    print("Leaf node insertion scenario triggered")
                    if not insert_leaf_at_terminal(current_node, insert_distance):
                        return
                else:
                    print(f"Checking insertion between previous node '{prev_node.name if prev_node else 'None'}' and current node '{current_node.name}' with distances {prev_dist} - {insert_distance}")
                    if not robust_insert_leaf_at_node(prev_node, prev_dist - insert_distance, current_node, prev_dist, toward_root):
                        return
                continue

            for child in current_node.children:
                if child not in visited_nodes:
                    print(f"Adding child node '{child.name}' to the queue")
                    queue.append((child, current_dist + child.dist, current_node, child.dist, current_path, False))

            if current_node.up and current_node.up not in visited_nodes:
                print(f"Adding parent node '{current_node.up.name}' to the queue")
                queue.append((current_node.up, current_dist + current_node.dist, current_node, current_node.dist, current_path, True))

    def validate_insertion_path(current_node, new_internal_node, previous_node, original_branch_distance):
        # Verifies if the insertion happened between the correct nodes
        print(f"Verifying insertion path...")
        distance_check = current_node.get_distance(new_internal_node) + new_internal_node.get_distance(previous_node)
        print(f"Verifying insertion path distance: {distance_check}, between '{previous_node.name}' and '{current_node.name}'")
        return abs(distance_check - original_branch_distance) < tolerance

    if dist <= target_node.dist:
        print(f"\nDirect insertion at target leaf '{target_leaf}' with distance {dist}")
        insert_leaf_at_terminal(target_node, dist)
    else:
        bfs(target_node, 0)

    def round_tree_distances(tree_node, decimals=8):
        for node in tree_node.traverse():
            node.dist = round(node.dist, decimals)

    round_tree_distances(tree)

    if insertion_points:
        print("\nFinal tree with all inserted leaves:")
        print(tree.write(format=1))
        print(tree)
    else:
        print("No valid insertion points were found based on the specified distance.")

# Example
newick = "(((A:1.587,(F:1.110,(M:1.343,R:1.369):0.846):0.487):1.981,D:0.356):2.121,(B:1.936,(C:0.915,Q:1.201):2.101):0.912);"
target_leaf = "Q"
new_leaf_base_name = "temp"
new_length = 0.279
dist = 3.0597060866386405

insert_leaf_from_target(newick, target_leaf, new_leaf_base_name, new_length, dist)
