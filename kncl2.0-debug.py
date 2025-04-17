# Debug version
"""
k-Nearest Common Leaves (k-NCL) Tree Completion Algorithm

This script implements the k-NCL algorithm for completing phylogenetic trees
that are defined on different but overlapping taxon sets. The implementation
uses the ete3 library to work with phylogenetic trees.
"""

from ete3 import Tree

######################################
# Helper function
######################################
def ensure_unique_internal_node_names(tree, prefix="INTERNAL"):
    """
    Ensures that every internal node in 'tree' has a unique name.
    If a node does not have a name or its name starts with the given prefix,
    it is renamed with a unique identifier.
    """
    counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf():
            if not node.name or node.name.startswith(prefix):
                node.name = f"{prefix}_{counter}"
                counter += 1

def check_leaves_in_tree(tree, leaf_names, tree_label="Tree"):
    """
    Check that each leaf name in 'leaf_names' exists in 'tree'. If some names
    are missing, print a warning.
    """
    missing = []
    for leaf in leaf_names:
        if not tree.search_nodes(name=leaf):
            missing.append(leaf)
    if missing:
        print(f"[WARNING] The following leaf names were not found in {tree_label}: {missing}")
    return len(missing)

#########################################
# Basic set operations for leaves
#########################################
def find_common_leaves(T1, T2):
    """
    Identify the set of common leaves between trees T1 and T2.
    """
    leaves1 = set(T1.get_leaf_names())
    leaves2 = set(T2.get_leaf_names())
    return leaves1.intersection(leaves2)

def find_distinct_leaves(T, common_leaves):
    """
    Return the set of leaves in tree T that are not in the set of common leaves.
    """
    all_leaves = set(T.get_leaf_names())
    return all_leaves - common_leaves

##########################################
# Finding maximal distinct-leaf subtrees
##########################################
def mark_all_distinct(node, distinct_leaf_set):
    """
    Mark node.allDistinct True if all descendant leaves belong to the distinct set,
    otherwise False. For leaves, check whether the leaf name is in distinct_leaf_set.
    """
    if node.is_leaf():
        node.add_feature("allDistinct", node.name in distinct_leaf_set)
    else:
        # Initially assume true, then check children.
        all_children_distinct = True
        for child in node.children:
            # Recursively mark child if not set.
            if not hasattr(child, "allDistinct"):
                mark_all_distinct(child, distinct_leaf_set)
            if not child.allDistinct:
                all_children_distinct = False
        node.add_feature("allDistinct", all_children_distinct)

def findSD(tree, distinct_leaves):
    """
    Identify maximal distinct-leaf subtrees in a phylogenetic tree.

    This function marks each node with the attribute `allDistinct`, which is True
    if all of the node’s descendant leaves belong to the given distinct leaf set (DL).
    It then returns a list of subtree roots where each subtree contains only distinct
    leaves and is maximal, i.e., no ancestor of the root is also entirely composed
    of distinct leaves.

    Parameters:
        tree (ete3.Tree): The phylogenetic tree to search.
        distinct_leaves (set of str): Set of distinct leaf names (DL).

    Returns:
        list of ete3.TreeNode: Roots of the maximal distinct-leaf subtrees (SD).
    """

    # Phase 1: Mark nodes with 'allDistinct' = True if all descendant leaves are in DL
    for node in tree.traverse("postorder"):
        mark_all_distinct(node, distinct_leaves)

    # Phase 2: Collect only the maximal distinct subtrees
    subtree_roots = []
    for node in tree.traverse("postorder"):
        if node.allDistinct:
            parent = node.up
            # Only include node if parent is not allDistinct (or no parent at all)
            if parent is None or not getattr(parent, "allDistinct", False):
                subtree_roots.append(node)

    return subtree_roots

#########################################
# Global and leaf-based adjustment rates
#########################################
def compute_global_adjustment_rate(T1, T2, common_leaves, debug_label="T1->T2"):
    """
    Compute the global adjustment rate r(T1, T2) as the ratio of the cumulative
    pairwise distances between all common leaves in T1 and T2.
    """
    check_leaves_in_tree(T1, common_leaves, f"T1 ({debug_label})")
    check_leaves_in_tree(T2, common_leaves, f"T2 ({debug_label})")

    cl_list = list(common_leaves)
    numerator = 0.0
    denominator = 0.0
    for i in range(len(cl_list)):
        for j in range(i + 1, len(cl_list)):
            l_i = cl_list[i]
            l_j = cl_list[j]
            d1 = T1.get_distance(l_i, l_j)
            d2 = T2.get_distance(l_i, l_j)
            numerator += d1
            denominator += d2

    if abs(denominator) < 1e-15:
        return 1.0
    return numerator / denominator

def compute_leaf_based_adjustment_rate(T_target, T_source, leaf_c, common_leaves):
    """
    Compute the leaf-based adjustment rate r^(l_c)(T_target, T_source) as the ratio:
        sum_{l_i in CL} d^(T_target)(l_c, l_i) / sum_{l_i in CL} d^(T_source)(l_c, l_i)
        """
    if not T_target.search_nodes(name=leaf_c) or not T_source.search_nodes(name=leaf_c):
        print(f"[WARNING] Leaf '{leaf_c}' not found in one of the trees during leaf-based rate calculation.")
        return 1.0

    numerator = 0.0
    denominator = 0.0
    for other_leaf in common_leaves:
        target_nodes = T_target.search_nodes(name=other_leaf)
        source_nodes = T_source.search_nodes(name=other_leaf)
        if not target_nodes or not source_nodes:
            print(f"[INFO] Skipping {other_leaf} in leaf-based rate: not found in both trees.")
            continue

        d1 = T_target.get_distance(leaf_c, target_nodes[0])
        d2 = T_source.get_distance(leaf_c, source_nodes[0])
        numerator += d1
        denominator += d2

    if abs(denominator) < 1e-15:
        return 1.0
    return numerator / denominator

##################################
# Subtree scaling
##################################
def scale_subtree(subtree_root, scale_factor):
    """
    Scale all branch lengths in subtree 'subtree_root' by 'scale_factor'.
    """
    for node in subtree_root.traverse():
        node.dist *= scale_factor

###################################
# k-nearest common leaves
###################################
def get_k_nearest_common_leaves(T, subtree_root, common_leaves, k):
    """
    Given tree T and a subtree with root at 'subtree_root', compute the k-nearest
    common leaves (by distance from the subtree root). If the k-th distance is tied,
    include all leaves with the same distance.
    """
    dist_list = []
    for l in common_leaves:
        if not T.search_nodes(name=l):
            continue
        d_val = T.get_distance(subtree_root, l)
        dist_list.append((l, d_val))
    dist_list.sort(key=lambda x: x[1])

    if not dist_list:
        return []

    result = []
    count = 0
    cutoff_dist = None
    for (leaf_name, d_val) in dist_list:
        if count < k:
            result.append(leaf_name)
            cutoff_dist = d_val
            count += 1
        else:
            if abs(d_val - cutoff_dist) < 1e-15:
                result.append(leaf_name)
            else:
                break

    return result

###################################################
# Minimizing the objective function for insertion
###################################################

def find_optimal_insertion_point(target_tree, subtree_root, ncl_list, leaf_based_d_p, target_distinct):
    """
    For each branch in target_tree (parameterized by x in [0,1)), this function computes
    an observed distance for each common leaf in ncl_list and evaluates the quadratic objective:
    
      OF(x) = sum_{l_c in ncl_list} (d(l_c, v(x)) - d_p(l_c, subtree_root))^2.
    
    Returns a tuple: (best_edge, best_x, best_obj)
    """
    best_edge = None
    best_x = 0.0
    best_obj = float('inf')
    MIN_TERMINAL_LENGTH = 1e-3

    # Ensure unique names for consistent matching.
    for node in target_tree.traverse("postorder"):
        if not node.name:
            node.name = f"AUTO_{id(node)}"
    
    for node in target_tree.traverse("postorder"):
        if node.up is None:
            continue  # Skip root
        
        # Exclude candidate edges whose descendant leaves intersect with target_distinct.
        descendant_leaves = set(node.get_leaf_names())
        if descendant_leaves.intersection(target_distinct):
            # Skip this candidate edge.
            continue

        # Also skip if the node or its parent has been used for a previous insertion.
        if getattr(node, "usedInsertion", False) or getattr(node.up, "usedInsertion", False):
            continue

        parent = node.up
        parent_candidates = target_tree.search_nodes(name=parent.name)
        if not parent_candidates:
            continue
        parent_in_target = parent_candidates[0]
        edge_length = node.dist
        if edge_length < 1e-15:
            continue

        A_vals = []
        sum_diff = 0.0
        for l_c in ncl_list:
            target_nodes = target_tree.search_nodes(name=l_c)
            if not target_nodes:
                print(f"[DEBUG] Common leaf '{l_c}' not in target tree; skipping.")
                continue
            dist_parent_lc = target_tree.get_distance(parent_in_target, target_nodes[0])
            T_val = leaf_based_d_p.get(l_c, 0.0)
            A_vals.append((dist_parent_lc, T_val))
            sum_diff += (T_val - dist_parent_lc)

        if not A_vals:
            continue

        # Compute candidate x_opt from the derivative zero condition.
        x_opt = sum_diff / (edge_length * len(A_vals))
        if node.is_leaf():
            if edge_length < MIN_TERMINAL_LENGTH:
                continue
            x_max = 1 - (MIN_TERMINAL_LENGTH / edge_length)
            x_clamped = max(0.0, min(x_opt, x_max))
        else:
            x_clamped = max(0.0, min(x_opt, 1.0))

        for x_test in [0.0, x_clamped]:
            of_val = 0.0
            for (A, T_val) in A_vals:
                observed = A + x_test * edge_length
                of_val += (observed - T_val) ** 2
            if of_val < best_obj:
                best_obj = of_val
                best_edge = (parent_in_target, node)
                best_x = x_test

        if best_edge is not None:
            print(f"[DEBUG] Candidate edge: Parent({parent_in_target.name}) -> Child({node.name}), "
                  f"edge_length={edge_length:.4f}, best_x={best_x:.4f}, OF={best_obj:.6f}")
    return best_edge, best_x, best_obj

def insert_subtree_at_point(target_tree, edge, x_opt, subtree_root):
    """
    Inserts 'subtree_root' into 'target_tree' at the candidate point along the edge
    defined by (parent, child) using the optimal x value. This version safely handles
    replacement of the original edge to avoid leaf duplication bugs.
    """
    parent, child = edge
    original_len = child.dist
    if abs(original_len) < 1e-15:
        parent.add_child(subtree_root)
        return

    dist_to_new = x_opt * original_len
    dist_to_child = (1 - x_opt) * original_len

    if x_opt < 1e-15:
        parent.add_child(subtree_root)
        return

    # Create a new internal node to split the branch
    new_node = Tree()
    new_node.dist = dist_to_new

    # Attempt to remove the child safely based on name match
    matching_children = [c for c in parent.children if c.name == child.name]
    if matching_children:
        for c in matching_children:
            parent.remove_child(c)
    else:
        print(f"[WARNING] Could not find exact child to remove: {child.name}")
        # Still add subtree directly to parent if removal fails
        parent.add_child(subtree_root)
        return

    # Proceed with edge splitting and insertion
    parent.add_child(new_node)

    # Reattach original child with adjusted branch length
    child.dist = dist_to_child
    new_node.add_child(child)

    # Attach the subtree
    connection_length = getattr(subtree_root, "original_root_branch_length", 0.0)
    subtree_root.dist = connection_length
    new_node.add_child(subtree_root)

    # Optional sanity check
    leaves = target_tree.get_leaf_names()
    dupes = [l for l in set(leaves) if leaves.count(l) > 1]
    if dupes:
        print(f"[DUPLICATE ALERT] Duplicated leaves after insertion: {dupes}")

def remove_internal_node_names(tree):
    """
    Remove names from all internal nodes in the tree to clean up the final output.
    """
    for node in tree.traverse():
        if not node.is_leaf():
            node.name = ""

########################################
# Main k-NCL function (with debug prints)
########################################
def kNCL(T1_original, T2_original, k=None):
    """
    Main function for k-Nearest Common Leaves tree completion.

    Steps:
      1) Make deep copies of the input trees.
      2) Assign unique names to all internal nodes.
      3) Identify common leaves between T1 and T2.
      4) Determine distinct leaves and extract maximal distinct-leaf subtrees.
      5) Compute global adjustment rates.
      6) Insert subtrees into each other's tree using the k-NCL logic.
    Returns:
        A tuple containing the completed trees (T1_plus, T2_plus).
    """
    print("[DEBUG] Starting k-NCL algorithm")

    T1 = T1_original.copy()
    T2 = T2_original.copy()

    ensure_unique_internal_node_names(T1, prefix="T1IN")
    ensure_unique_internal_node_names(T2, prefix="T2IN")

    # 1) Identify common leaves. Require at least two common leaves.
    CL = find_common_leaves(T1, T2)
    print(f"[DEBUG] Common leaves ({len(CL)}): {sorted(CL)}")
    if len(CL) < 2:
        raise ValueError("Need at least two common leaves to proceed.")

    # Set default k as the number of common leaves if not specified
    if k is None:
        k = len(CL)

    # 2) Identify distinct leaves in each tree.
    DL1 = find_distinct_leaves(T1, CL)
    DL2 = find_distinct_leaves(T2, CL)
    print(f"[DEBUG] Distinct leaves in T1: {sorted(DL1)}")
    print(f"[DEBUG] Distinct leaves in T2: {sorted(DL2)}")

    # 3) Extract maximal distinct-leaf subtrees.
    SD1 = findSD(T1, DL1)
    SD2 = findSD(T2, DL2)
    print(f"[DEBUG] Identified {len(SD1)} distinct-leaf subtrees in T1")
    for idx, st in enumerate(SD1, 1):
        print(f"  T1 subtree {idx} root: {st.name}, leaves: {sorted(st.get_leaf_names())}")
    print(f"[DEBUG] Identified {len(SD2)} distinct-leaf subtrees in T2")
    for idx, st in enumerate(SD2, 1):
        print(f"  T2 subtree {idx} root: {st.name}, leaves: {sorted(st.get_leaf_names())}")

    # 4) Compute global adjustment rates.
    r12 = compute_global_adjustment_rate(T1, T2, CL, debug_label="T1->T2")
    r21 = compute_global_adjustment_rate(T2, T1, CL, debug_label="T2->T1")
    #print(f"[DEBUG] Global adjustment rate r(T1->T2): {r12}")
    #print(f"[DEBUG] Global adjustment rate r(T2->T1): {r21}")
    
    def insert_subtrees_into_target(target_tree, source_tree, subtree_roots, global_rate, label, CL):
        for root_s in subtree_roots:
            original_root_branch = root_s.dist if root_s.up is not None else 0.0
            subtree_copy = root_s.copy()
            subtree_copy.add_feature("original_root_branch_length", original_root_branch * global_rate)
            scale_subtree(subtree_copy, global_rate)
            print(f"\n[DEBUG] Inserting subtree (root: {root_s.name}) into {label}")
            #print(f"  - Original root branch length: {original_root_branch}")
            #print(f"  - Scaled root branch length: {subtree_copy.original_root_branch_length}")

            max_k = len(CL)
            k_current = min(k, max_k)
            unique_insertion_found = False

            while k_current <= max_k:
                ncl = get_k_nearest_common_leaves(source_tree, root_s, CL, k_current)
                #print(f"  - k-Nearest common leaves (k={k_current}): {ncl}")

                if not ncl:
                    print(f"  - No common leaves found for k={k_current}, skipping insertion.")
                    break

                leaf_based_d_p = {}
                for l_c in ncl:
                    r_lc = compute_leaf_based_adjustment_rate(target_tree, source_tree, l_c, CL)
                    dist_source = source_tree.get_distance(l_c, root_s) if source_tree.search_nodes(name=l_c) else 0.0
                    d_p_val = dist_source * r_lc
                    leaf_based_d_p[l_c] = d_p_val
                    #print(f"    - Leaf: {l_c}, Dist(source): {dist_source:.4f}, r^(l_c): {r_lc:.4f}, d_p: {d_p_val:.4f}")

                target_distinct = find_distinct_leaves(target_tree, CL)
                edge, x_opt, best_obj = find_optimal_insertion_point(target_tree, root_s, ncl, leaf_based_d_p, target_distinct)


                #edge, x_opt, best_obj = find_optimal_insertion_point(target_tree, root_s, ncl, leaf_based_d_p)

                if edge is not None:
                    unique_insertion_found = True
                    break

                print(f"  [DEBUG] No unique insertion found for k={k_current}, retrying with k={k_current + 1}")
                k_current += 1

            if not unique_insertion_found:
                print("  - Insertion ambiguous or failed after retries, attaching at root.")
                target_tree.add_child(subtree_copy)
            else:
                #print(f"  - Inserting at edge: ({edge[0].name}, {edge[1].name}), x_opt: {x_opt:.4f}, ObjVal: {best_obj:.6f}")
                insert_subtree_at_point(target_tree, edge, x_opt, subtree_copy)

            # Mark all nodes in the inserted subtree so they are not reused
            for node in subtree_copy.traverse():
                node.add_feature("isInserted", True)

            # Ensure all internal nodes have unique names again
            ensure_unique_internal_node_names(target_tree, prefix="TINS")

        return target_tree
    
    # 5) Insert subtrees from T2 into T1
    print("\n[DEBUG] Inserting T2 subtrees into T1")
    T1_completed = insert_subtrees_into_target(T1, T2, SD2, r12, label="T1", CL=CL)

    # 6) Insert subtrees from T1 into T2
    print("\n[DEBUG] Inserting T1 subtrees into T2")
    T2_completed = insert_subtrees_into_target(T2, T1, SD1, r21, label="T2", CL=CL)

    # 7) Clean final trees
    remove_internal_node_names(T1_completed)
    remove_internal_node_names(T2_completed)

    # 8) Final validation: Confirm both trees have the same leaf set and match the union of L(T1) and L(T2)
    final_leaves_T1 = set(T1_completed.get_leaf_names())
    final_leaves_T2 = set(T2_completed.get_leaf_names())
    original_union = set(T1_original.get_leaf_names()) | set(T2_original.get_leaf_names())

    print("\n[VALIDATION] Final leaf count check:")
    print(f"  - L(T1) = {len(T1_original.get_leaf_names())}")
    print(f"  - L(T2) = {len(T2_original.get_leaf_names())}")
    print(f"  - L(T1) ∪ L(T2) = {len(original_union)}")
    print(f"  - Leaves in completed T1 = {len(final_leaves_T1)}")
    print(f"  - Leaves in completed T2 = {len(final_leaves_T2)}")

    if final_leaves_T1 == final_leaves_T2 == original_union:
        print("[SUCCESS] Completed trees have consistent and complete leaf sets.")
    else:
        print("[ERROR] Leaf sets in completed trees do not match L(T1) ∪ L(T2)!")
        missing_T1 = original_union - final_leaves_T1
        missing_T2 = original_union - final_leaves_T2
        extra_T1 = final_leaves_T1 - original_union
        extra_T2 = final_leaves_T2 - original_union
        if missing_T1:
            print(f"  - Missing in T1+: {sorted(missing_T1)}")
        if missing_T2:
            print(f"  - Missing in T2+: {sorted(missing_T2)}")
        if extra_T1:
            print(f"  - Extra in T1+: {sorted(extra_T1)}")
        if extra_T2:
            print(f"  - Extra in T2+: {sorted(extra_T2)}")

    print("[DEBUG] Tree completion finished\n")
    return T1_completed, T2_completed
