## Cluster construction methodology

For each of the four groups (Amphibians, Birds, Mammals, and Sharks), we start with five base phylogenetic trees whose leaf sets overlap by 10–90%. Each base tree defines a cluster. From every base tree, four additional trees are simulated by introducing biological processes (e.g., horizontal gene transfer) using the [AsymmeTree library](https://github.com/david-schaller/AsymmeTree) ([Schaller et al., 2022](https://doi.org/10.3390/software1030013)). In total, this yields 25 trees per species group, organized as five consecutive clusters of five trees.

For every unique tree pair in each cluster dataset, we run two completion strategies: k-NCL and RF(+). After completion, we compute RF distances for trees completed with RF(+), and for k-NCL–completed trees we compute both RF(k-NCL) and BSD(k-NCL) distances. We then assemble pairwise distance matrices for all three metrics (RF(+), RF(k-NCL), BSD(k-NCL)) and visualize them with heatmaps and boxplots to assess how clearly each completion method recovers the intended cluster structure. 

This folder includes tree cluster datasets for four species groups, precomputed distance matrices, and a script that draws the heatmaps and boxplots from those matrices.
