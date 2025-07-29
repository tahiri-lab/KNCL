# *k*-NCL
***k*-Nearest Common Leaves: overlapping phylogenetic tree completion**

> The results of this work will be submitted.
>
> Koshkarov, A. & Tahiri, N. (2025). *k*-Nearest Common Leaves algorithm for phylogenetic tree completion (in progress).

#### Description :bookmark_tabs:
`kncl.py` is a Python script designed to complete 2 phylogenetic trees defined on different but mutually overlapping sets of taxa on the union of their taxa sets. The script processes trees in Newick format by pairwise comparing them. It outputs the completed versions of these trees.

>:pushpin: **Current Status**
>
>Please note that `kncl.py` is an active work in progress. Development is ongoing to improve its efficiency and extend its functionality to cover a wider range of use cases and tree types (rooted and unrooted).

  - Biologically meaningful datasets of partially overlapping phylogenetic trees with branch lengths are available [here](https://github.com/tahiri-lab/KNCL/tree/main/data/).
  - Biological datasets used in the evaluation part of the *k*-NCL algorithm are located [here](https://github.com/tahiri-lab/KNCL/tree/main/data/evaluation_datasets).
  - The data construction pipeline is available [here](https://github.com/tahiri-lab/KNCL/tree/main/data/data-pipeline).

#### Requirements :clipboard:
- Python 3.x
- `ete3` Python package (with its required dependencies)

#### Installation :wrench:
Ensure Python 3 and `ete3` are installed. You can install `ete3` via pip if it's not already installed:
```bash
pip install ete3
```

#### Usage :bulb:
Run the script using the command:
```bash
python kncl.py -i <input.newick> <k> -o <output.txt>
```
- `<input.newick>`: File containing two or more trees in Newick format, each tree on a separate line.
- `<k>`: Integer value of *k* for the *k*-nearest common leaves algorithm (*k* must be between 2 and the number of common leaves).
- `<output.txt>`: File where the output will be saved.

#### Example :bookmark:
Given an input file `example.newick` with the following content:
```
((((Hylobates:4.711,Nomascus:4.711):7.489,Pongo:12.100):5.946,Macaca:18.050):9.590,(Ateles:10.699,Alouatta:10.699):17.037);
((((Hylobates:5.420,Nomascus:5.420):6.377,((Pan:4.945,Gorilla:4.945):4.227,Pongo:9.172):2.624):7.520,Macaca:19.317):8.377,Ateles:27.693);
```
Run the script as:
```bash
python kncl.py -i example.newick 2 -o output.txt
```

#### Output :book:
The `output.txt` file will contain:
- Completed versions of each tree pair.

#### Notes :pencil:
- Ensure the input Newick file is correctly formatted.

The description of the *k*-Nearest Common Leaves algorithm can be found [here](https://github.com/tahiri-lab/KNCL/blob/main/algorithm.md).

# 📧 Contact
Please email us at: <Nadia.Tahiri@USherbrooke.ca> for any questions or feedback.
