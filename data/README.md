## Datasets

This folder contains datasets of phylogenetic trees defined on different but overlapping sets of taxa.

### Simulated data

Simulated phylogenetic tree data were generated using the GPTree generator (sourses: [1](https://github.com/tahiri-lab/GPTree), [2](https://github.com/tahiri-lab/GPTree/tree/GPTreeCluster)), which is specifically designed to simulate phylogenetic trees with varying numbers of leaves and different levels of overlap. The dataset was produced by configuring the generator with parameters that included generating a total of 500 trees, setting the number of leaves to range between 20 and 30, and adjusting the overlap levels to 0.3, 0.4, 0.5, 0.6, and 0.7 (100 trees for each level of overlap).

### Biological data

The methodology for obtaining biological data of phylogenetic trees with different but overlapping taxa is as follows.
1. The biological data utilized in this study was obtained from [vertlife.org](https://vertlife.org/phylosubsets/), which offers a straightforward method for acquiring tree distributions with specified subsets of taxa. The tool initially prunes a comprehensive dataset to a smaller subset and then samples trees from the selected pseudoposterior distribution.
2. The following four groups, representing four distinct datasets, have been selected for analysis: **amphibians**, **birds**, **mammals**, and **sharks**.
3. The number of species included in each group varies. In particular, there are 7239 species of amphibians, 9993 species of birds, 5911 species of mammals, and 1192 species of sharks (see [all_species_lists.xlsx](https://github.com/tahiri-lab/KNCL/blob/main/data/all_species_lists.xlsx)).

To be continued...
