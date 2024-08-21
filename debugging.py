#Debugging of the updated version of the k-NCL algorithm

from ete3 import Tree

# Helper functions
def precompute_descendants(node, distinct_leaves):
    if node.is_leaf():
        node.add_feature("descendants_distinct", node.name in distinct_leaves)
    else:
        all_distinct = True
        for child in node.children:
            if not child.descendants_distinct:
                all_distinct = False
                break
        node.add_feature("descendants_distinct", all_distinct)

def findSD(tree, distinct_leaves):
    for node in tree.traverse("postorder"):
        precompute_descendants(node, distinct_leaves)

    subtree_roots = set()
    visited_leaves = set()

    for leaf_name in distinct_leaves:
        if leaf_name in visited_leaves:
            continue

        leaf_node = tree & leaf_name
        current_node = leaf_node
        subtree_root = None

        while current_node:
            if current_node.descendants_distinct:
                subtree_root = current_node
            current_node = current_node.up

        if subtree_root and subtree_root not in subtree_roots:
            subtree_roots.add(subtree_root)
            visited_leaves.update(get_leaves(subtree_root))

    return subtree_roots

def get_leaves(node):
    return set(leaf.name for leaf in node)

def get_subtree_newick_with_branch_lengths(node):
    def recursive_newick(node):
        if node.is_leaf():
            return f"{node.name}:{node.dist}"
        children_newick = [recursive_newick(child) for child in node.children]
        return f"({','.join(children_newick)}):{node.dist}"

    return recursive_newick(node)

def InsertTempLeaves(newick, target_leaf, new_leaf_base_name, new_length, dist, tolerance=1e-10):
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
        if abs(insert_distance) < tolerance:
            #new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            parent_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(parent_node)
            visited_nodes.add(parent_node)
            visited_nodes.add(new_leaf_name)
        else:
            excess_length = parent_node.dist - insert_distance
            if excess_length < 0:
                return False

            if parent_node.up:
                new_internal_node = parent_node.up.add_child(dist=excess_length)
                parent_node.detach()
                new_internal_node.add_child(parent_node, dist=insert_distance)
                #new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
                new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
                new_internal_node.add_child(name=new_leaf_name, dist=new_length)
                insertion_points.append(new_internal_node)
                visited_nodes.add(new_internal_node)
                visited_nodes.add(new_internal_node.children[1])  # Newly added leaf node
            else:
                # Handling the case when the parent node is the root
                new_internal_node = tree.add_child(dist=excess_length)
                parent_node.detach()
                new_internal_node.add_child(parent_node, dist=insert_distance)
                #new_leaf_name = f"{new_leaf_base_name}{len(insertion_points) + 1}"
                new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
                new_internal_node.add_child(name=new_leaf_name, dist=new_length)
                insertion_points.append(new_internal_node)
                visited_nodes.add(new_internal_node)
                visited_nodes.add(new_internal_node.children[1])  # Newly added leaf node
            return True

    def bfs_case_1(node, accumulated_distance):
        queue = [(node, accumulated_distance, None, 0, [])]
        while queue:
            current_node, current_dist, prev_node, prev_dist, path = queue.pop(0)
            if current_node in visited_nodes:
                continue
            visited_nodes.add(current_node)
            current_path = path + [current_node.name]

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
        insert_leaf_at_node(target_node, dist, target_node)
    else:
        # Determine which BFS method to use based on tree topology
        if tree.get_distance(tree.get_tree_root(), target_node) < tree.get_farthest_leaf()[1] / 2:
            bfs_case_1(target_node, 0)
        else:
            bfs_case_2(target_node, 0)

    return tree

def find_farthest_leaf(tree, start, temporary_leaves):
    max_distance = 0
    farthest_leaf = start
    for leaf in temporary_leaves:
        if leaf is not start:
            distance = tree.get_distance(start, leaf)
            if distance > max_distance:
                max_distance = distance
                farthest_leaf = leaf
    return farthest_leaf, max_distance

def find_path(leaf1, leaf2):
    path1 = []
    node = leaf1
    while node:
        path1.append(node)
        node = node.up

    path2 = []
    node = leaf2
    while node:
        path2.append(node)
        node = node.up

    lca = None
    for n1 in path1:
        if n1 in path2:
            lca = n1
            break

    path = []
    branch_lengths = []
    for node in path1:
        path.append(node)
        if node == lca:
            break

    lca_index = path2.index(lca)
    reversed_path2 = list(reversed(path2[:lca_index]))
    path.extend(reversed_path2)

    for i in range(1, len(path)):
        branch_lengths.append(path[i - 1].get_distance(path[i]))

    path_names = [n.name for n in path]
    return path, branch_lengths

def compute_midpoint(tree, temporary_leaves):
    start = next(iter(temporary_leaves))
    leaf1, dist1 = find_farthest_leaf(tree, start, temporary_leaves)
    leaf2, dist2 = find_farthest_leaf(tree, leaf1, temporary_leaves)
    path, branch_lengths = find_path(leaf1, leaf2)
    total_distance = dist2
    half_distance = round(total_distance / 2, 8)

    cumulative_distance = 0
    prev_node = None
    for i, node in enumerate(path):
        if i > 0:
            cumulative_distance += round(branch_lengths[i - 1], 8)
        if cumulative_distance >= half_distance:
            excess = round(cumulative_distance - half_distance, 8)
            prev_node = path[i - 1]
            return prev_node, node, excess, half_distance, branch_lengths[i - 1]
        elif cumulative_distance == half_distance:
            prev_node = path[i - 1]
            return prev_node, node, 0, half_distance, branch_lengths[i - 1]

def insert_midpoint_and_new_leaf(tree, prev_node, curr_node, excess, new_leaf_name, branch_length, original_dist):
    if excess == 0:
        new_leaf = Tree(name=new_leaf_name)
        new_leaf.dist = branch_length
        curr_node.add_child(new_leaf)
        return tree

    if curr_node in prev_node.get_ancestors():
        distance_to_midpoint = round(excess, 8)
        distance_from_midpoint_to_leaf = round(original_dist - excess, 8)
    else:
        distance_to_midpoint = round(original_dist - excess, 8)
        distance_from_midpoint_to_leaf = round(excess, 8)

    if distance_to_midpoint < 0 or distance_from_midpoint_to_leaf < 0:
        raise ValueError("Negative distance encountered. Check the calculation logic.")

    new_node = Tree(name="midpoint")
    new_node.dist = distance_to_midpoint

    if curr_node in prev_node.get_ancestors():
        parent = prev_node.up
        child = prev_node
    else:
        parent = prev_node
        child = curr_node

    parent.remove_child(child)
    parent.add_child(new_node)
    new_node.add_child(child)
    child.dist = distance_from_midpoint_to_leaf

    new_leaf = Tree(name=new_leaf_name)
    new_leaf.dist = branch_length
    new_node.add_child(new_leaf)

    return tree

def remove_temporary_leaves(tree, temporary_leaves):
    for leaf in temporary_leaves:
        if leaf.up:
            parent = leaf.up
            parent.remove_child(leaf)

def kNCL(T1, T2, k):
    CL = set(T1.get_leaf_names()) & set(T2.get_leaf_names())
    print(f"Common Leaves (CL): {CL}")

    def adjust_rate(T1, T2):
        sum_T1 = sum(T1.get_distance(l1, l2) for i, l1 in enumerate(CL) for l2 in list(CL)[i + 1:])
        sum_T2 = sum(T2.get_distance(l1, l2) for i, l1 in enumerate(CL) for l2 in list(CL)[i + 1:])
        return sum_T1 / sum_T2 if sum_T2 else 1

    r12 = adjust_rate(T1, T2)
    r21 = adjust_rate(T2, T1)
    print(f"Adjustment rates: r12 = {r12}, r21 = {r21}")

    SD1 = findSD(T1, set(T1.get_leaf_names()) - CL)
    SD2 = findSD(T2, set(T2.get_leaf_names()) - CL)
    print(f"Subtrees in T1 (SD1): {[get_subtree_newick_with_branch_lengths(n) for n in SD1]}")
    print(f"Subtrees in T2 (SD2): {[get_subtree_newick_with_branch_lengths(n) for n in SD2]}")

    def process_tree(T1, T2, SD1, r12, k):
        for a in SD1:
            if not a.name:
                a.name = "subtree_" + str(len(SD1))  # Assign a name to unnamed subtrees
            br_a = a.dist * r12
            print(f"Processing subtree {a.name} with branch length {br_a}")
            NCL = sorted(CL, key=lambda l: T1.get_distance(a, T1 & l))[:k]
            print(f"Nearest Common Leaves for {a.name}: {NCL}")
            TL = set()
            for lc in NCL:
                print(f"Checking for leaf {lc} in T2")
                if T2.search_nodes(name=lc):
                    lc_node = T2 & lc
                else:
                    raise ValueError(f"Common leaf '{lc}' not found in T2.")
                print(f"Checking for leaf {lc} in T1")
                if T1.search_nodes(name=lc):
                    lc_node_T1 = T1 & lc
                else:
                    raise ValueError(f"Common leaf '{lc}' not found in T1.")
                rc = sum(T2.get_distance(lc_node, l) for l in CL) / sum(T1.get_distance(lc_node_T1, l) for l in CL)
                dp = (T1.get_distance(a, lc_node_T1) - a.dist) * rc
                print(f"Inserting temporary leaves for {a.name} from leaf {lc} at distance {dp}")
                T2 = InsertTempLeaves(T2.write(format=1), lc, "temp", br_a, dp)
            print(f"Temporary leaves for {a.name}: {[leaf.name for leaf in T2.get_leaves() if 'temp' in leaf.name]}")
            print("Tree after inserting temporary leaves:")
            print(T2.write(format=1))
            temporary_leaves = [T2 & leaf_name for leaf_name in T2.get_leaf_names() if 'temp' in leaf_name]
            prev_node, curr_node, excess, _, original_dist = compute_midpoint(T2, temporary_leaves)
            print(f"Inserting midpoint and new leaf for {a.name} between {prev_node.name} and {curr_node.name}")
            T2 = insert_midpoint_and_new_leaf(T2, prev_node, curr_node, excess, a.name, br_a, original_dist)
            remove_temporary_leaves(T2, temporary_leaves)
            print(f"Inserted midpoint and new leaf for {a.name}")
            print(f"Tree after removing temporary leaves:")
            print(T2.write(format=1))
        return T2

    T1_completed = process_tree(T1, T2, SD1, r12, k)
    print("Intermediate T1 completed tree state:")
    print(T1_completed.write(format=1))

    T2_completed = process_tree(T2, T1, SD2, r21, k)
    print("Intermediate T2 completed tree state:")
    print(T2_completed.write(format=1))

    return T1_completed, T2_completed

# Example
newick1 = "(((L4:1.1058,L6:0.7225):0.6678,(L10:0.5582,(L1:0.8540,L5:0.6621):1.1539):1.1164):1.8280,(L3:0.6057,L8:1.0215):0.9939,(L7:1.1467,(2:1.0002,9:1.3349):1.5883):0.6401);"
newick2 = "(((L1:0.5,L2:1.0):0.5,L3:1.5):0.5,(L4:1.0,(L5:1.5,L6:2.0):1.0):0.5);"

T1 = Tree(newick1, format=1)
T2 = Tree(newick2, format=1)

k = 3

T1_completed, T2_completed = kNCL(T1, T2, k)

print("Completed T1:")
print(T1_completed.write(format=1))
print("Completed T2:")
print(T2_completed.write(format=1))
