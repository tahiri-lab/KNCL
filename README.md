# K-NCL
**k-Nearest Common Leaves: phylogenetic tree completion and distance calculation**

> The results of this work will be submitted.
>
> Koshkarov, A. & Tahiri, N. (2024). Novel algorithm for comparing phylogenetic trees with different but overlapping taxa, (in progress).

#### Description :bookmark_tabs:
`bsd_distance.py` is a Python script designed to complete 2 phylogenetic trees defined on different but mutually overlapping sets of taxa on the union of their taxa sets and compute the Branch Score Distance (BSD(+)) between the completed phylogenetic trees. The BSD is a measure of the dissimilarity between two trees based on their branch lengths. The script processes trees in Newick format, comparing the first tree in the input file with every other tree. It outputs the completed versions of these trees, their BSD(+) distance, the pruned versions of the trees, and the BSD(-) distance.

>:pushpin: **Current Status**
>
>Please note that the current version of `bsd_distance.py` is under development. It is not a final version and requires further improvements and testing to ensure robustness and accuracy in various use cases.

#### Requirements :clipboard:
- Python 3.x
- `ete3` Python package

#### Installation :wrench:
Ensure Python 3 and `ete3` are installed. You can install `ete3` via pip if it's not already installed:
```bash
pip install ete3
```

#### Usage :bulb:
Run the script using the command:
```bash
python3 bsd_distance.py -i <input.newick> <k> -o <output.txt>
```
- `<input.newick>`: File containing two or more trees in Newick format, each tree on a separate line.
- `<k>`: Integer value of k for the k-nearest common leaves algorithm.
- `<output.txt>`: File where the output will be saved.

#### Example :bookmark:
Given an input file `example.newick` with the following content:
```
(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);
(E:0.1,F:0.2,(C:0.3,D:0.4):0.5);
```
Run the script as:
```bash
python3 bsd_distance.py -i example.newick 2 -o output.txt
```

#### Output :book:
The `output.txt` file will contain:
- Completed versions of each tree pair.
- BSD(+) distances between the completed trees.
- Pruned versions of the trees on their common taxa.
- BSD(-) distances for the pruned trees.

#### Notes :pencil:
- Ensure the input Newick file is correctly formatted.
- The script handles trees with no common leaves by joining them under a new root.
- The BSD(+) calculations assume the trees are correctly parsed and all necessary functions are implemented.

The description of the k-Nearest Common Leaves algorithm can be found [here](https://github.com/tahiri-lab/KNCL/blob/main/algorithm.md).
