## Datasets

This folder contains datasets of phylogenetic trees defined on different but overlapping sets of taxa.

### Simulated data :computer:

Simulated phylogenetic tree data were generated using the GPTree generator (sources: [GPTree](https://github.com/tahiri-lab/GPTree) [[1]](#ref1), [GPTree Cluster](https://github.com/tahiri-lab/GPTree/tree/GPTreeCluster) [[2]](#ref2)), which is specifically designed to simulate phylogenetic trees with varying numbers of leaves and different levels of overlap. The dataset was produced by configuring the generator with parameters that included generating a total of 500 trees, setting the number of leaves to range between 20 and 30, and adjusting the overlap levels to 0.3, 0.4, 0.5, 0.6, and 0.7 (100 trees for each level of overlap).

The simulated data consists of the following two files:

1. [**`simulated_data_gptree.txt`**](https://github.com/tahiri-lab/KNCL/blob/main/data/simulated_data_gptree.txt):
   - Contains 500 phylogenetic trees simulated by [GPTree](https://github.com/tahiri-lab/GPTree) [1].
   - Includes 100 trees for each overlap level: 30%, 40%, 50%, 60%, and 70%.

2. [**`simulated_data_gptree_cluster.txt`**](https://github.com/tahiri-lab/KNCL/blob/main/data/simulated_data_gptree_cluster.txt):
   - Contains 500 phylogenetic trees simulated by [GPTreeCluster](https://github.com/tahiri-lab/GPTree/tree/GPTreeCluster) [2].
   - Includes 5 clusters of 20 trees for each overlap level: 30%, 40%, 50%, 60%, and 70%.


### Biological data :deciduous_tree:

The methodology for obtaining biological data of phylogenetic trees with different but overlapping taxa is as follows.
1. The biological data utilized in this study was obtained from [vertlife.org](https://vertlife.org/phylosubsets/) [[3]](#ref3), which offers a straightforward method for acquiring tree distributions with specified subsets of taxa. The tool initially prunes a comprehensive dataset to a smaller subset and then samples trees from the selected pseudoposterior distribution.
2. The following four groups, representing four distinct datasets, have been selected for analysis: :frog: **amphibians** [[4]](#ref4), :eagle: **birds** [[5]](#ref5), :monkey: **mammals** [[3]](#ref3), and :shark: **sharks** [[6]](#ref6).
3. The number of species included in each group varies. In particular, there are 7239 species of amphibians, 9993 species of birds, 5911 species of mammals, and 1192 species of sharks (see [all_species_lists.xlsx](https://github.com/tahiri-lab/KNCL/blob/main/data/all_species_lists.xlsx)).
4. Each species group comprises subgroups representing species families. For example, species such as Acris blanchardi, Acris crepitans, and Acris gryllus are assumed to represent the Acris subgroup, which is identified by the first word in the species names. In each subgroup, a single species is randomly selected.
5. The resulting aggregated dataset for each species subgroup is subsequently employed in the creation of a dataset of overlapping trees. The final dataset for each group consists of a number of subsets that overlap from 10% to 90% (see the following picture).

![Levels of overlap among subsets for 4 groups of species](https://github.com/tahiri-lab/KNCL/blob/main/data/images/overlaps_subsets.png "Levels of overlap among subsets for 4 groups of species")

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

3. :monkey: [**`mammals_trees.txt`**](https://github.com/tahiri-lab/KNCL/blob/main/data/mammals_trees.txt):
   - Contains 800 phylogenetic trees (Newick) of *Mammals* with overlap levels ranging from 10% to 100%.
   - Total number of unique species: 104.
   - Average level of overlap: 58.85%.
   - Number of unique pairs of trees: 319600.
   - Number of unique pairs of trees with 100% overlap: 31600 or 9.89%.

4. :shark: [**`sharks_trees.txt`**](https://github.com/tahiri-lab/KNCL/blob/main/data/sharks_trees.txt):
   - Contains 900 phylogenetic trees (Newick) of *Sharks* with overlap levels ranging from 10% to 100%.
   - Total number of unique species: 69.
   - Average level of overlap: 59.21%.
   - Number of unique pairs of trees: 404550.
   - Number of unique pairs of trees with 100% overlap: 40050 or 9.9%.

### References

1. <a id="ref1"></a> Koshkarov, A., & Tahiri, N. (2023). GPTree: Generator of Phylogenetic Trees with Overlapping and Biological Events for Supertree Inference. In *BIOINFORMATICS* (pp. 212-219). [https://doi.org/10.5220/0011697100003414](https://doi.org/10.5220/0011697100003414)

2. <a id="ref2"></a> Koshkarov, A., & Tahiri, N. (2023). GPTree Cluster: phylogenetic tree cluster generator in the context of supertree inference. *Bioinformatics Advances*, 3(1). [https://doi.org/10.1093/bioadv/vbad023](https://doi.org/10.1093/bioadv/vbad023)

3. <a id="ref3"></a> Upham, N. S., J. A. Esselstyn, and W. Jetz. 2019. Inferring the mammal tree: species-level sets of phylogenies for questions in ecology, evolution, and conservation. *PLOS Biology*. [https://doi.org/10.1371/journal.pbio.3000494](https://doi.org/10.1371/journal.pbio.3000494)

4. <a id="ref4"></a> Jetz, W., and R. A. Pyron. 2018. The interplay of past diversification and evolutionary isolation with present imperilment across the amphibian tree of life. *Nature Ecology & Evolution*, 1. [https://www.nature.com/articles/s41559-018-0515-5](https://www.nature.com/articles/s41559-018-0515-5)

5. <a id="ref5"></a> Jetz, W., G. H. Thomas, J. B. Joy, K. Hartmann, and A. O. Mooers. 2012. The global diversity of birds in space and time. *Nature*, 491:444–448. [http://www.nature.com/nature/journal/v491/n7424/abs/nature11631.html](http://www.nature.com/nature/journal/v491/n7424/abs/nature11631.html)

6. <a id="ref6"></a> Stein, R. W., Mull, C. G., Kuhn, T. S., Aschliman, N. C., Davidson, L. N. K., Joy, J. B., Smith, G. J., Dulvy, N. K., & Mooers, A. O. 2018. Global priorities for conserving the evolutionary history of sharks, rays, and chimaeras. *Nature Ecology & Evolution*, 2:288–298. [http://dx.doi.org/10.1038/s41559-017-0448-4](http://dx.doi.org/10.1038/s41559-017-0448-4)
