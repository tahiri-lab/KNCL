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
    return node.write(format=1)

def InsertTempLeaves(tree, target_leaf, new_leaf_base_name, new_length, dist, tolerance=1e-10):
    # Operate directly on 'tree'
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
            insertion_points.append(new_leaf_name)
            visited_nodes.add(new_internal_node)
        else:
            new_internal_node = parent.add_child(dist=excess_length)
            new_internal_node.add_child(previous_node, dist=insert_distance)
            new_leaf_name = f"{target_leaf}_{new_leaf_base_name}{len(insertion_points) + 1}"
            new_internal_node.add_child(name=new_leaf_name, dist=new_length)
            insertion_points.append(new_leaf_name)
            visited_nodes.add(new_internal_node)

        return True

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
            insertion_points.append(new_leaf_name)
            visited_nodes.add(new_internal_node)
        else:
            return False

        return True

    def bfs(node, accumulated_distance):
        queue = [(node, accumulated_distance, None, 0, [], False)]
        while queue:
            current_node, current_dist, prev_node, prev_dist, path, toward_root = queue.pop(0)
            if current_node in visited_nodes or 'temp' in current_node.name:
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

    if dist <= target_node.dist:
        insert_leaf_at_terminal(target_node, dist)
    else:
        bfs(target_node, 0)

    return insertion_points  # Return names instead of node objects

def find_farthest_leaf(tree, start, temporary_leaves):
    max_distance = 0
    farthest_leaf = start
    for leaf_name in temporary_leaves:
        if leaf_name != start.name:
            leaf = tree & leaf_name
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

    return path, branch_lengths

def compute_midpoint(tree, temporary_leaves):
    start_name = next(iter(temporary_leaves))
    start = tree & start_name
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

def insert_midpoint_and_new_subtree(tree, prev_node, curr_node, excess, subtree, branch_length, original_dist):
    if excess == 0:
        # Attach the subtree directly to current node curr_node
        new_subtree = subtree.copy()
        new_subtree.dist = branch_length
        curr_node.add_child(new_subtree)
        return tree

    if curr_node in prev_node.get_ancestors():
        distance_to_midpoint = round(excess, 8)
        distance_from_midpoint_to_leaf = round(original_dist - excess, 8)
    else:
        distance_to_midpoint = round(original_dist - excess, 8)
        distance_from_midpoint_to_leaf = round(excess, 8)

    if distance_to_midpoint < 0 or distance_from_midpoint_to_leaf < 0:
        raise ValueError("Negative distance encountered. Check the calculation logic.")

    new_node = Tree()
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

    # Now add the subtree
    new_subtree = subtree.copy()
    new_subtree.dist = branch_length
    new_node.add_child(new_subtree)

    return tree

def remove_temporary_leaves(tree, temporary_leaves):
    def collapse_single_child_nodes(node):
        while not node.is_leaf() and len(node.children) == 1:
            child = node.children[0]
            child.dist += node.dist
            parent = node.up
            if parent:
                parent.remove_child(node)
                parent.add_child(child)
                node = parent
            else:
                # Node is root
                child.up = None
                tree = child  # Update tree reference
                node = child
    for leaf_name in temporary_leaves:
        leaf_nodes = tree.search_nodes(name=leaf_name)
        if leaf_nodes:
            leaf = leaf_nodes[0]
            parent = leaf.up
            if parent:
                parent.remove_child(leaf)
                collapse_single_child_nodes(parent)
            else:
                # Leaf is root
                tree = None
    return tree  # Return the possibly updated tree

def clear_internal_node_names(tree):
    for node in tree.traverse():
        if not node.is_leaf():
            node.name = ''

def kNCL(T1, T2, k):
    CL = set(T1.get_leaf_names()) & set(T2.get_leaf_names())
    if len(CL) < 3:
        raise ValueError("The input trees must have at least 3 common leaves.")
    if k < 2 or k > len(CL):
        raise ValueError("The value of k must be between 2 and the number of common leaves.")

    print(f"Common Leaves (CL): {CL}")

    def adjust_rate(T1, T2):
        sum_T1 = sum(T1.get_distance(l1, l2) for i, l1 in enumerate(CL) for l2 in list(CL)[i + 1:])
        sum_T2 = sum(T2.get_distance(l1, l2) for i, l1 in enumerate(CL) for l2 in list(CL)[i + 1:])
        return sum_T1 / sum_T2 if sum_T2 else 1

    r12 = adjust_rate(T1, T2)  # Adjusting from T2 to T1
    r21 = adjust_rate(T2, T1)  # Adjusting from T1 to T2
    print(f"Adjustment rates: r12 = {r12}, r21 = {r21}")

    SD1 = findSD(T1, set(T1.get_leaf_names()) - CL)
    SD2 = findSD(T2, set(T2.get_leaf_names()) - CL)
    print(f"Subtrees in T1 (SD1): {[get_subtree_newick_with_branch_lengths(n) for n in SD1]}")
    print(f"Subtrees in T2 (SD2): {[get_subtree_newick_with_branch_lengths(n) for n in SD2]}")

    def process_tree(target_tree, source_tree, subtrees_to_insert, rate, k):
        for a in subtrees_to_insert:
            if not a.name:
                a.name = "subtree_" + str(len(subtrees_to_insert))  # Assign a name to unnamed subtrees

            # Make a copy of the subtree to avoid modifying the original
            adjusted_subtree = a.copy()

            # Adjust all branch lengths in the copied subtree
            for node in adjusted_subtree.traverse():
                node.dist *= rate

            print(f"Processing subtree {a.name} with adjusted branch lengths")
            NCL = sorted(CL, key=lambda l: source_tree.get_distance(a, source_tree & l))[:k]
            print(f"Nearest Common Leaves for {a.name}: {NCL}")
            TL = set()
            for lc in NCL:
                print(f"Checking for leaf {lc} in target_tree")
                if target_tree.search_nodes(name=lc):
                    lc_node = target_tree & lc
                else:
                    raise ValueError(f"Common leaf '{lc}' not found in target_tree.")
                print(f"Checking for leaf {lc} in source_tree")
                lc_node_source = source_tree & lc
                rc = sum(target_tree.get_distance(lc_node, l) for l in CL) / sum(source_tree.get_distance(lc_node_source, l) for l in CL)
                dp = (source_tree.get_distance(a, lc_node_source) - a.dist) * rc
                print(f"Inserting temporary leaves for {a.name} from leaf {lc} at distance {dp}")
                temp_leaves = InsertTempLeaves(target_tree, lc, "temp", adjusted_subtree.dist, dp)
                print(f"Temporary leaves after insertion for {lc}: {temp_leaves}")
                TL.update(temp_leaves)

            print(f"Temporary leaves for {a.name}: {TL}")
            print("Tree after inserting temporary leaves:")
            print(target_tree.write(format=1))

            if not TL:
                print(f"No temporary leaves were inserted for {a.name}")
                continue

            if len(TL) == 1:
                single_leaf_name = next(iter(TL))
                single_leaf = target_tree & single_leaf_name
                parent = single_leaf.up
                branch_length = single_leaf.dist

                # Remove the temporary leaf
                parent.remove_child(single_leaf)

                # Attach the adjusted subtree to the parent
                adjusted_subtree.dist = branch_length
                parent.add_child(adjusted_subtree)

                print(f"Only one temporary leaf, inserted subtree {a.name}")
            else:
                try:
                    prev_node, curr_node, excess, _, original_dist = compute_midpoint(target_tree, TL)
                    print(f"Inserting midpoint and new subtree for {a.name} between {prev_node.name} and {curr_node.name}")
                    print(f"Midpoint insertion details - Excess: {excess}, Original Dist: {original_dist}")
                    target_tree = insert_midpoint_and_new_subtree(target_tree, prev_node, curr_node, excess, adjusted_subtree, adjusted_subtree.dist, original_dist)
                except Exception as e:
                    print("Error encountered during midpoint insertion:")
                    print(f"Nodes involved: {TL}")
                    print(f"Error: {str(e)}")

            target_tree = remove_temporary_leaves(target_tree, TL)
            print(f"Inserted midpoint and new subtree for {a.name}")
            print(f"Tree after removing temporary leaves:")
            print(target_tree.write(format=1))
        return target_tree

    # Process trees with the corrected adjustment rates
    T1_completed = process_tree(T1, T2, SD2, r12, k)
    print("Intermediate T1 completed tree state:")
    print(T1_completed.write(format=1))

    T2_completed = process_tree(T2, T1, SD1, r21, k)
    print("Intermediate T2 completed tree state:")
    print(T2_completed.write(format=1))

    clear_internal_node_names(T1_completed)
    clear_internal_node_names(T2_completed)

    return T1_completed, T2_completed

# Test example
newick1 = "((A:0.597,B:0.139):0.735,((C:0.171,E:0.069):0.218,(Q:0.138,D:0.077):0.343):0.609);"
newick2 = "(((A:1.587,(F:1.110,(M:1.343,R:1.369):0.846):0.487):1.981,D:0.356):2.121,(B:1.936,(C:0.915,Q:1.201):2.101):0.912);"

T1 = Tree(newick1, format=1)
T2 = Tree(newick2, format=1)

k = 3

T1_completed, T2_completed = kNCL(T1, T2, k)

print("\nCompleted T1:")
print(T1_completed.write(format=1))
print("Completed T2:")
print(T2_completed.write(format=1))
