# K-NCL
**k-Nearest Common Leaves: phylogenetic tree completion algorithm**

> The results of this work were submitted to [RECOMB-CG 2024](http://recomb-cg.org/), the 21st Annual Satellite Conference of RECOMB on Comparative Genomics: Koshkarov, A. & Tahiri, N. (2024). Novel algorithm for comparing phylogenetic trees with different but overlapping taxa. In RECOMB International Workshop on Comparative Genomics.

We have 2 phylogenetic trees in Newick format. These trees have different but overlapping taxa (leaves). The tree completion process consists of making both trees defined on the same taxa (this is the union of the leaf sets of both trees) by adding distinct (non-common) leaves from one tree to the other. This process involves the following five steps.

**Step 1.** For both trees find their common leaves, distinct (non-common) leaves and the maximal distict-leaf subtrees. A subtree $S$ of a tree $T$ is called a maximal distinct-leaf subtree (MDS) if and only if all leaves in the subtree $S$ belong to the set $DL(T)$ and there is no other subtree of $T$ that includes all the leaves from $DL(T)$ and has $S$ as its proper subtree. In other words, a maximal distinct-leaf subtree is a subtree of $T$ that includes only leaves from $DL(T)$ and cannot be further extended while still satisfying the condition of including only leaves from $DL(T)$. If a set of distinct leaves $DL(T)$ is available, then the elements from this set can be used as starting points. For each unvisited distinct leaf, we can perform a breadth-first traversal starting from that leaf to construct a Maximal Distinct-leaf Subtree (MDS).

For a tree $T$, the set containing all its maximal distinct-leaf subtrees and remaining distinct leaves is denoted as $SD(T)$.

**Step 2.** Calculate the following:
1. distances (without repetitions) between the common leaves for both trees (based on the definition 2). The distance between any two leaves $l_1$ and $l_2$ of the phylogenetic tree T is the sum of branch lengths on the unique path from $l_1$ to $l_2$.

2. the branch adjustment rates for both trees. The branch adjustment rate is the ratio of the sums of pairwise (without repetitions) distances between common leaves in one tree to the other. The rates are defined by the following formulas:

$r(T_1, T_2) = \frac{\sum\limits_{i=1}^{N_{CL}-1} \sum\limits_{j=i+1}^{N_{CL}}d^{(T_1)}(l_i,l_j)}{\sum\limits_{i=1}^{N_{CL}-1} \sum\limits_{j=i+1}^{N_{CL}}d^{(T_2)}(l_i,l_j)}$,

$r(T_2, T_1) = \frac{1}{r(T_1, T_2)}$, where $N_{CL}$ is the number of common leaves $CL(T_1,T_2)$.

3. the leaf-based adjustment rates for each common leaf. For any common leaf $l_c \in CL(T_1,T_2)$, the leaf-based adjustment rates for trees $T_1$ and $T_2$ are defined by the following formulas:

$r^{(l_c)}(T_1, T_2) = \frac{\sum\limits_{i=1}^{N_{CL}-1}d^{(T_1)}(l_c,l_i)}{\sum\limits_{i=1}^{N_{CL}-1}d^{(T_2)}(l_c,l_i)}$,

$r^{(l_c)}(T_2, T_1) = \frac{\sum\limits_{i=1}^{N_{CL}-1}d^{(T_2)}(l_c,l_i)}{\sum\limits_{i=1}^{N_{CL}-1}d^{(T_1)}(l_c,l_i)}$,

where $r^{(l_c)}(T_1, T_2)$ is the l_c-based adjustment rate for tree $T_1$ related to tree $T_2$, $r^{(l_c)}(T_2, T_1)$ is the $l_c$-based adjustment rate for tree $T_2$ related to tree $T_1$, $l_c, l_i \in CL(T_1,T_2)$, $l_i \neq l_c$. Each leaf-based adjustment rate is calculated based on one common leaf, relative to the other common leaves in the considered trees.

**Step 3.** Calculate the new root branch lengths for each maximal distinct-leaf subtrees (the branch connecting the root of a maximal distinct-leaf subtree $S$ to its lowest ancestor node in the tree $T$ is called the root branch of that subtree) and the new terminal branch for the remaining distinct leaves as the current branch lengths multiplied by the corresponding branch adjustment rate, $r(T_2, T_1)$ for selected branches in tree 1 and $r(T_1, T_2)$ for selected branches using the following equations:

$\forall a \in SD(T_1), br^{(T_2)}(a) = br^{(T_1)}(a) \cdot r(T_2, T_1)$,

$\forall b \in SD(T_2), br^{(T_1)}(b) = br^{(T_2)}(b) \cdot r(T_1, T_2)$.

All branch lengths within maximal distinct-leaf subtrees should be adjusted according to the appropriate rate when they are transferred from one tree to another. It is important to note that the initial branch lengths in the trees should be preserved without any changes, and we change the length of the branches on their copies. In the following steps we will insert those leaves and maximal distinct-leaf subtrees (if any) with their calculated branch lengths from one tree to another.

**Step 4.** For each element inside $SD(T_1)$ for tree $T_1$ and inside $SD(T_2)$ for tree $T_2$ do the following:

1. Select k nearest common leaves $l_k$ (the user needs to input this value) sorted in ascending order by distance to the considered element from SD and calculate the cutback distances between each selected common leaf and that element. The cutback distance between any two leaves $l_1$ and $l_2$ of the phylogenetic tree $T$ (denoted as $dc^{(T)}(l_1,l_2)$ ) is determined by the sum of the branch lengths between $l_1$ and the parent node of leaf $l_2$. Alternatively, the cutback distance $dc^{(T)}(l_1,l_2)$ in the phylogenetic tree $T$ is calculated as the distance between leaves $l_1$ and $l_2$ minus the length of the terminal branch leading to leaf $l_2$:

$dc^{(T)}(l_1,l_2) = d^{(T)}(l_1,l_2) - br^{(T)}(l_2)$,
where $br^{(T)}(l_2)$ is the terminal branch leading to leaf $l_2$ in tree $T$.

3. Multiply these distances by the corresponding leaf-based adjustment rates (the two following formulas) to obtain the distances $d_p$ in order to find possible positions for adding new leaves related to the same common leaves in another tree:
    
$\forall a \in SD(T_1), \forall l_k \in CL(T_1,T_2), d_p^{(T_2)}(l_k,a) = dc^{(T_1)}(l_k,a) \cdot r^{(l_k)}(T_2, T_1)$,

$\forall b \in SD(T_2), \forall l_k \in CL(T_1,T_2), d_p^{(T_1)}(l_k,b) = dc^{(T_2)}(l_k,b) \cdot r^{(l_k)}(T_1, T_2)$.

3. Add temporary nodes in all possible positions (only among the branches that were in the tree initially, not including new branches and leaves) at appropriate calculated distances ($d_p$) from the same common leaves in the second tree. That is, the selected nearest common leaves $l_k$ in the second tree are starting points for new temporary nodes, at a computed distance from which they are inserted.

4. When all such temporary nodes are added, find the midpoint among temporary nodes and planting the considered distinct leaf (or maximal distinct-leaf subtree) with its new branch length at this position. In the case where no temporary nodes can be found for a distinct leaf or a maximal distinct-leaf subtree within the tree topology, the considered distinct leaf (or maximal distinct-leaf subtree) with its new branch length should be planted at the midpoint of $k$ selected common leaves.

**Step 5.** As a result, we obtain 2 phylogenetic trees $T_1^{\uplus}$ and $T_2^{\uplus}$ defined on the same set of taxa (on the union of the sets of leaves of the initial trees $L(T_1) \cup L(T_2)$ ). The following formula is used to calculate the distance between these trees:

$BSD(+)(T_1^{\uplus},T_2^{\uplus}) = \sqrt{\sum\limits_{i=1}^{N_{\cup}-1} \sum\limits_{j=i+1}^{N_{\cup}}(d^{(T_1^{\uplus})}(l_i,l_j)-d^{(T_2^{\uplus})}(l_i,l_j))^2}$,

where $d^{(T_1^{\uplus})}(l_i,l_j)$ is the distance between leaves $l_i$ and $l_j$ in tree $T_1^{\uplus}$, $d^{(T_2^{\uplus})}(l_i,l_j)$ is the distance between leaves $l_i$ and $l_j$ in tree $T_2^{\uplus}$, $l_i, l_j \in L(T_1) \cup L(T_2)$, $N_{\cup}$ is the size of the set $L(T_1) \cup L(T_2)$.
