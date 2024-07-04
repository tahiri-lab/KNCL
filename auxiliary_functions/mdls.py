# Finding maximal distinct-leaf subtrees

from ete3 import Tree

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

def find_maximal_distinct_subtrees(tree, distinct_leaves):
    # Precompute whether all descendant leaves are distinct for each node
    for node in tree.traverse("postorder"):
        precompute_descendants(node, distinct_leaves)

    subtree_roots = set()
    visited_leaves = set()

    for leaf_name in distinct_leaves:
        if leaf_name in visited_leaves:
            continue

        leaf_node = tree&leaf_name
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

def get_subtree_newick_without_branch_lengths(node):
    def recursive_newick(node):
        if node.is_leaf():
            return node.name
        children_newick = [recursive_newick(child) for child in node.children]
        return f"({','.join(children_newick)})"

    return recursive_newick(node)

# Example
if __name__ == "__main__":
    # Create a sample phylogenetic tree using ete3
    newick = "(((L4:1.1058,L6:0.7225):0.6678,(L10:0.5582,(L1:0.8540,L5:0.6621):1.1539):1.1164):1.8280,(L3:0.6057,L8:1.0215):0.9939,(L7:1.1467,(2:1.0002,9:1.3349):1.5883):0.6401);"
    tree = Tree(newick, format=1)

    distinct_leaves = {"L1", "L3", "L4", "L5", "L7", "L8", "L10"}

    result = find_maximal_distinct_subtrees(tree, distinct_leaves)
    for subtree_root in result:
        subtree_newick = get_subtree_newick_without_branch_lengths(subtree_root)
        print(subtree_newick)
