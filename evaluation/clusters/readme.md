## Cluster construction methodology

For each of the four groups (Amphibians, Birds, Mammals, and Sharks), we start with five base phylogenetic trees whose leaf sets overlap by 10–90%. Each base tree defines a cluster. From every base tree, four additional trees are simulated by introducing biological processes (e.g., horizontal gene transfer) using the [AsymmeTree library](https://github.com/david-schaller/AsymmeTree) ([Schaller et al., 2022](https://doi.org/10.3390/software1030013)). In total, this yields 25 trees per species group, organized as five consecutive clusters of five trees.
