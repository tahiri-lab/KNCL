# K-NCL
**k-Nearest Common Leaves: phylogenetic tree completion algorithm**

We have 2 phylogenetic trees in Newick format. These trees have different but overlapping taxa (leaves). The tree completion process consists of making both trees defined on the same taxa (this is the union of the leaf sets of both trees) by adding distinct (non-common) leaves from one tree to the other. 

The $k$-Nearest Common Leaves algorithm can be summarized at a high level as follows. The algorithm begins by identifying the common leaves, distinct leaves, and maximal distinct-leaf subtrees in both phylogenetic trees. Next, it calculates the branch adjustment rates to ensure that evolutionary distances are preserved when integrating distinct leaves and subtrees. Temporary leaves are then placed at specific positions within the tree to facilitate this integration. The algorithm computes midpoints among the temporary leaves to determine the insertion points for the distinct leaves. Finally, the distinct leaves are inserted into the tree at these midpoints, with branch lengths adjusted proportionally to maintain the evolutionary relationships and distances. The pseudocode of the algorithm is shown below.

![k-Nearest Common Leaves Algorithm](https://github.com/tahiri-lab/KNCL/blob/main/img/kncl.png "k-NCL")

Auxiliary functions can be found [here](https://github.com/tahiri-lab/KNCL/tree/main/auxiliary_functions).

As a result, we obtain 2 phylogenetic trees $T_1^{\uplus}$ and $T_2^{\uplus}$ defined on the same set of taxa (on the union of the sets of leaves of the initial trees $L(T_1) \cup L(T_2)$ ). The following formula is used to calculate the distance between these trees:

$BSD(+)(T_1^{\uplus},T_2^{\uplus}) = \sqrt{\sum\limits_{i=1}^{N_{\cup}-1} \sum\limits_{j=i+1}^{N_{\cup}}(d^{(T_1^{\uplus})}(l_i,l_j)-d^{(T_2^{\uplus})}(l_i,l_j))^2}$,

where $d^{(T_1^{\uplus})}(l_i,l_j)$ is the distance between leaves $l_i$ and $l_j$ in tree $T_1^{\uplus}$, $d^{(T_2^{\uplus})}(l_i,l_j)$ is the distance between leaves $l_i$ and $l_j$ in tree $T_2^{\uplus}$, $l_i, l_j \in L(T_1) \cup L(T_2)$, $N_{\cup}$ is the size of the set $L(T_1) \cup L(T_2)$.
