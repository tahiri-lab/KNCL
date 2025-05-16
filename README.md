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
- `ete3` Python package

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
(((((((Callimic:0.023313,Ateles:0.010045)59:0.003307,(((Pan:0.001001,Homo:0.003222)98:0.007021,Nomascus:0.019337)65:0.006297,Macaca:0.022545)75:0.003800)100:0.056141,Tarsius:0.070541)40:0.004811,(Otolemur:0.080291,Lemur:0.073501)67:0.014141)54:0.014589,Tupaia:0.110178)85:0.046160,Cynoceph:0.040415)100:0.356615,Rattus:0.048351,Mus:0.036439);
(((((Hylobate:0.006938,(Pongo:0.007054,(Pan:0.002692,Gorilla:0.003234)100:0.002698)96:0.001954)100:0.008626,Macaca:0.020688)100:0.010416,(Ateles:0.009814,Alouatta:0.013133)100:0.024863)100:0.030640,(Tarsius:0.097437,(Otolemur:0.075033,Lemur:0.046192)100:0.019542)65:0.007328)100:0.185167,Rattus:0.048223,Mus:0.063981);
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
