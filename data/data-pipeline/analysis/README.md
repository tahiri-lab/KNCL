# Analysis

This folder contains datasets and scripts used for evaluating the performance and functionality of the proposed method for assembling phylogenetic datasets with overlapping taxa, as described in the article "New method for assembling biological datasets of phylogenetic trees with overlapping taxa." The analyses include supertree construction results and runtime evaluations.

This includes:

- **`input_multisets/`**

  Contains 100 input sets of phylogenetic trees used for the supertree reconstruction demonstration. Each file (`multiset_X.txt`) consists of 30 trees with overlapping taxa, generated using the proposed pipeline.

- **`supertrees/`**

  This folder includes output supertrees produced by five different methods:
  - `supertrees_sfit.txt` — Split Fit algorithm
  - `supertrees_dfit.txt` — Most Similar Supertree
  - `supertrees_nj.txt` — Average NJ
  - `supertrees_mrplus.txt` — Majority-Rule
  - `supertrees_scs.txt` — Spectral Clustering

- **`correctness_validation.ipynb`**

This script performs correctness validation for phylogenetic tree datasets generated using the data pipeline.

- **`runtime_evaluation_for_data_pipeline.ipynb`**

This script includes results related to measuring the runtime performance of the dataset assembly pipeline under various configurations (e.g., number of trees, taxon overlap, etc.).

- **`supertree_validation_for_data_pipeline.ipynb`** 

  Python script that processes the above supertrees and input sets to calculate:
  - Success rate of reconstruction;
  - Average number of taxa in the output trees;
  - Average Robinson-Foulds (RF) distance between each supertree and its corresponding input trees.

