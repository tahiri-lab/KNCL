# K-NCL
**k-Nearest Common Leaves: phylogenetic tree completion algorithm**

We have 2 phylogenetic trees in Newick format. These trees have different but overlapping taxa (leaves). The tree completion process consists of making both trees defined on the same taxa (this is the union of the leaf sets of both trees) by adding distinct (non-common) leaves from one tree to the other. 

The $k$-Nearest Common Leaves algorithm can be summarized at a high level as follows. The algorithm begins by identifying the common leaves, distinct leaves, and maximal distinct-leaf subtrees in both phylogenetic trees. Next, it calculates the branch adjustment rates to ensure that evolutionary distances are preserved when integrating distinct leaves and subtrees. Temporary leaves are then placed at specific positions within the tree to facilitate this integration. The algorithm computes midpoints among the temporary leaves to determine the insertion points for the distinct leaves. Finally, the distinct leaves are inserted into the tree at these midpoints, with branch lengths adjusted proportionally to maintain the evolutionary relationships and distances. The pseudocode of the algorithm is shown below.

![k-Nearest Common Leaves Algorithm](https://github.com/tahiri-lab/KNCL/blob/main/img/kncl.png "k-NCL")

Auxiliary functions can be found [here](https://github.com/tahiri-lab/KNCL/tree/main/auxiliary_functions).

**Necessary formulas:**

The branch adjustment rate $r(T_2, T_1)$ is defined as the ratio of the sum of pairwise distances among the common leaves in $T_1$ to the sum of the pairwise distances among the same leaves in $T_2$:

$r(T_2, T_1) = \frac{\sum\limits_{l_i, l_j \in CL(T_1, T_2), i < j} d^{(T_2)}(l_i, l_j)}{\sum\limits_{l_i, l_j \in CL(T_1, T_2), i < j} d^{(T_1)}(l_i, l_j)}$.

$br^{(T_2^{\uplus})}(a) = br^{(T_1)}(a) \cdot r(T_2, T_1)$, where $a$ represents an element from $SD(T_1)$.

$r^{(l_c)}(T_2, T_1) = \frac{\sum\limits_{i=1}^{n_{CL}}d^{(T_2)}(l_c,l_i)}{\sum\limits_{i=1}^{n_{CL}}d^{(T_1)}(l_c,l_i)}$, where $r^{(l_c)}(T_2, T_1)$ is the adjustment rate specific to leaf $l_c$, and $l_c, l_i \in CL(T_1,T_2)$.

The temporary position distance for an element $a \in SD(T_1)$, which is to be newly added and is associated with a specified common leaf $l_c$ in tree $T_2$, is then utilized to identify possible locations for temporary leaves within the tree topology in the tree completion process. This distance is computed as follows.

$d_p^{(T_2^{\uplus})}(l_c,a) = \left(d^{(T_1)}(l_c,a) - br^{(T_1)}(a) \right) \cdot r^{(l_c)}(T_2, T_1)$.

As a result, we obtain 2 phylogenetic trees $T_1^{\uplus}$ and $T_2^{\uplus}$ defined on the same set of taxa (on the union of the sets of leaves of the initial trees $L(T_1) \cup L(T_2)$ ). The following formula is used to calculate the distance between these trees:

$BSD(+)(T_1^{\uplus},T_2^{\uplus}) = \sqrt{\sum\limits_{i=1}^{N_{\cup}-1} \sum\limits_{j=i+1}^{N_{\cup}}(d^{(T_1^{\uplus})}(l_i,l_j)-d^{(T_2^{\uplus})}(l_i,l_j))^2}$,

where $d^{(T_1^{\uplus})}(l_i,l_j)$ is the distance between leaves $l_i$ and $l_j$ in tree $T_1^{\uplus}$, $d^{(T_2^{\uplus})}(l_i,l_j)$ is the distance between leaves $l_i$ and $l_j$ in tree $T_2^{\uplus}$, $l_i, l_j \in L(T_1) \cup L(T_2)$, $N_{\cup}$ is the size of the set $L(T_1) \cup L(T_2)$.
