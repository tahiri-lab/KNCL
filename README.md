# K-NCL
**k-Nearest Common Leaves: phylogenetic tree completion and distance calculation**

> The results of this work will be submitted.
>
> Koshkarov, A. & Tahiri, N. (2024). *k*-Nearest Common Leaves algorithm for phylogenetic tree completion (in progress).

#### Description :bookmark_tabs:
`kncl.py` is a Python script designed to complete 2 phylogenetic trees defined on different but mutually overlapping sets of taxa on the union of their taxa sets and compute the Branch Score Distance (BSD(+)) between the completed phylogenetic trees. The BSD is a measure of the dissimilarity between two trees based on their branch lengths (see more details [here](https://www.mdpi.com/2073-8994/16/7/790)). The script processes trees in Newick format, comparing the first tree in the input file with every other tree. It outputs the completed versions of these trees, their BSD(+) distance, and the BSD(-) distance.

>:pushpin: **Current Status**
>
>Please note that the current version of `kncl.py` is under development. It is not a final version and requires further improvements and testing to ensure robustness and accuracy in various use cases.

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
python3 kncl.py -i <input.newick> <k> -o <output.txt>
```
- `<input.newick>`: File containing two or more trees in Newick format, each tree on a separate line.
- `<k>`: Integer value of *k* for the *k*-nearest common leaves algorithm (*k* must be between 2 and the number of common leaves).
- `<output.txt>`: File where the output will be saved.

#### Example :bookmark:
Given an input file `example.newick` with the following content:
```
(((((((Callimic:0.023313,Ateles:0.010045)59:0.003307,(((Pan:0.001001,Homo:0.003222)98:0.007021,Nomascus:0.019337)65:0.006297,Macaca:0.022545)75:0.003800)100:0.056141,Tarsius:0.070541)40:0.004811,(Otolemur:0.080291,Lemur:0.073501)67:0.014141)54:0.014589,Tupaia:0.110178)85:0.046160,Cynoceph:0.040415)100:0.356615,Rattus:0.048351,Mus:0.036439);
(((((Hylobate:0.006938,(Pongo:0.007054,(Pan:0.002692,Gorilla:0.003234)100:0.002698)96:0.001954)100:0.008626,Macaca:0.020688)100:0.010416,(Ateles:0.009814,Alouatta:0.013133)100:0.024863)100:0.030640,(Tarsius:0.097437,(Otolemur:0.075033,Lemur:0.046192)100:0.019542)65:0.007328)100:0.185167,Rattus:0.048223,Mus:0.063981);
```
Run the script as:
```bash
python3 kncl.py -i example.newick 2 -o output.txt
```

#### Output :book:
The `output.txt` file will contain:
- BSD(+) distances between the completed trees.
- BSD(-) distances for the pruned trees.
- Completed versions of each tree pair.

#### Notes :pencil:
- Ensure the input Newick file is correctly formatted.

The description of the *k*-Nearest Common Leaves algorithm can be found [here](https://github.com/tahiri-lab/KNCL/blob/main/algorithm.md).

Debugging of the updated version of the *k*-NCL algorithm is in progress.

# 📧 Contact
Please email us at: <Nadia.Tahiri@USherbrooke.ca> for any questions or feedback.
