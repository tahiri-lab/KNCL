## Datasets

This folder contains datasets of phylogenetic trees defined on different but overlapping sets of taxa.

### Simulated data :computer:

Simulated phylogenetic tree data were generated using the GPTree generator (sourses: [1](https://github.com/tahiri-lab/GPTree), [2](https://github.com/tahiri-lab/GPTree/tree/GPTreeCluster)), which is specifically designed to simulate phylogenetic trees with varying numbers of leaves and different levels of overlap. The dataset was produced by configuring the generator with parameters that included generating a total of 500 trees, setting the number of leaves to range between 20 and 30, and adjusting the overlap levels to 0.3, 0.4, 0.5, 0.6, and 0.7 (100 trees for each level of overlap).

The simulated data consists of the following two files:

1. [**`simulated_data_gptree.txt`**](https://github.com/tahiri-lab/KNCL/blob/main/data/simulated_data_gptree.txt):
   - Contains 500 phylogenetic trees simulated by [GPTree](https://github.com/tahiri-lab/GPTree).
   - Includes 100 trees for each overlap level: 30%, 40%, 50%, 60%, and 70%.

2. [**`simulated_data_gptree_cluster.txt`**](https://github.com/tahiri-lab/KNCL/blob/main/data/simulated_data_gptree_cluster.txt):
   - Contains 500 phylogenetic trees simulated by [GPTreeCluster](https://github.com/tahiri-lab/GPTree/tree/GPTreeCluster).
   - Includes 5 clusters of 20 trees for each overlap level: 30%, 40%, 50%, 60%, and 70%.


### Biological data :deciduous_tree:

The methodology for obtaining biological data of phylogenetic trees with different but overlapping taxa is as follows.
1. The biological data utilized in this study was obtained from [vertlife.org](https://vertlife.org/phylosubsets/), which offers a straightforward method for acquiring tree distributions with specified subsets of taxa. The tool initially prunes a comprehensive dataset to a smaller subset and then samples trees from the selected pseudoposterior distribution.
2. The following four groups, representing four distinct datasets, have been selected for analysis: **amphibians** :frog:, **birds** :eagle:, **mammals** :monkey:, and **sharks** :shark:.
3. The number of species included in each group varies. In particular, there are 7239 species of amphibians, 9993 species of birds, 5911 species of mammals, and 1192 species of sharks (see [all_species_lists.xlsx](https://github.com/tahiri-lab/KNCL/blob/main/data/all_species_lists.xlsx)).
4. Each species group comprises subgroups representing species families. For example, species such as Acris blanchardi, Acris crepitans, and Acris gryllus are assumed to represent the Acris subgroup, which is identified by the first word in the species names. In each subgroup, a single species is randomly selected.
5. The resulting aggregated dataset for each species subgroup is subsequently employed in the creation of a dataset of overlapping trees. The final dataset for each group consists of a number of subsets that overlap from 10% to 90%.

The biological data consists of the following four files:

1. :frog: [**`amphibians_trees.txt`**](https://github.com/tahiri-lab/KNCL/blob/main/data/amphibians_trees.txt):
   - Contains 700 phylogenetic trees (Newick) of *Amphibians* with overlap levels ranging from 10% to 100%.
   - Total number of unique species: 137.
   - Average level of overlap: 59.52%.
   - Number of unique pairs of trees: 244650.
   - Number of unique pairs of trees with 100% overlap: 24150 or 9.87%.
  
2. :eagle: [**`birds_trees.txt`**](https://github.com/tahiri-lab/KNCL/blob/main/data/birds_trees.txt):
   - Contains 600 phylogenetic trees (Newick) of *Birds* with overlap levels ranging from 10% to 100%.
   - Total number of unique species: 173.
   - Average level of overlap: 58.68%.
   - Number of unique pairs of trees: 179700.
   - Number of unique pairs of trees with 100% overlap: 17700 or 9.85%.

To be continued...
