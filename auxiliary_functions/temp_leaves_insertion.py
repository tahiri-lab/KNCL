# Inserting temporary leaves into a phylogenetic tree

import io
from Bio import Phylo
from Bio.Phylo.Newick import Clade

def insert_temp_leaves(tree, leaf_distances):
    def traverse_and_insert(clade, target_leaf, branch_length, distance, current_distance=0):
        new_clades = []
        for sub_clade in clade.clades:
            new_current_distance = current_distance + sub_clade.branch_length
            if sub_clade.is_terminal() and sub_clade.name == target_leaf:
                if distance < sub_clade.branch_length:
                    insertion_length = distance - current_distance
                    remaining_length = sub_clade.branch_length - insertion_length

                    new_internal_clade = Clade(branch_length=insertion_length)
                    temp_leaf = Clade(branch_length=branch_length, name=f"{target_leaf}_1")
                    new_internal_clade.clades.append(temp_leaf)
                    new_internal_clade.clades.append(Clade(branch_length=remaining_length, clades=[sub_clade]))

                    clade.clades.remove(sub_clade)
                    clade.clades.append(new_internal_clade)
                else:
                    temp_leaf = Clade(branch_length=branch_length, name=f"{target_leaf}_1")
                    clade.clades.append(temp_leaf)
                return True  # Stop recursion once the target leaf is found
            else:
                if traverse_and_insert(sub_clade, target_leaf, branch_length, distance, new_current_distance):
                    return True
        return False

    for leaf_info in leaf_distances:
        leaf, branch_length, distance = leaf_info.split(':')
        branch_length = float(branch_length)
        distance = float(distance)
        
        traverse_and_insert(tree.root, leaf, branch_length, distance)

    return tree

# Testing
newick_str = "((A:0.5,B:0.3):0.7,(C:0.6,D:0.4):0.8);"
tree = Phylo.read(io.StringIO(newick_str), "newick")

leaf_distances = ['A:0.6:1.22', 'B:0.4:0.2'] #branch lengths for temporary leaves and target distances from selected common leaves
modified_tree = insert_temp_leaves(tree, leaf_distances)

# Print the modified tree
Phylo.write(modified_tree, sys.stdout, "newick")
