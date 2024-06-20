## Phylogenetic Tree Branch Score Distance (BSD) Calculation

This script calculates the Branch Score Distance (BSD) between two phylogenetic trees provided in Newick format. Depending on the overlap of the leaf sets of the two trees, it calculates either the BSD or the BSD(-) value. The script handles three scenarios:

1. **Identical Leaf Sets**: Calculates the BSD.
2. **Overlapping but Different Leaf Sets**: Prunes the trees to their common leaves and calculates the BSD(-).
3. **No Common Leaves**: Informs the user that BSD cannot be computed.

### Usage

To use this script, run it from the command line with the required arguments:

```bash
python bsd.py -i input_file.txt -o output_file.txt
```

### Arguments

- `-i` or `--input`: Path to the input file containing two phylogenetic trees in Newick format (one per line).
- `-o` or `--output`: Path to the output file where the results will be written.


### Output

The output will be written to the specified output file. The script will provide the following results based on the leaf set overlap:

- **BSD:** When both trees have identical leaf sets.
- **BSD(-):** When the trees have overlapping but not identical leaf sets.
- An informative message when there are no common leaves.

### Additional details

1. **Parsing and Initial Analysis**: The script parses the Newick formatted trees and identifies the common leaves.
2. **Distance Calculations**: Calculates pairwise distances between the leaves of each tree.
3. **BSD Calculation**: 
   - If the trees have the same set of leaves, it calculates the BSD.
   - If the trees have overlapping but different sets of leaves, it prunes the trees to their common leaves and calculates the BSD(-).
   - If there are no common leaves, it prints a message indicating that the BSD cannot be computed.

### Example Output

For the example input provided, the output in `results.txt` might look like:

```
BSD: 0.8452
```

or

```
BSD(-): 1.0345
```

or

```
The BSD distance between input trees cannot be computed because these trees have no common leaves.
```

### Requirements

- Python 3.x
- `ete3` library

### Installation

- This script utilizes the `ete3` library for parsing and manipulating phylogenetic trees.

Install the required library using pip:

```bash
pip install ete3
```

### License

This project is licensed under the MIT License.
