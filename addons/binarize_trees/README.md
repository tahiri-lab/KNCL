# Binarize Newick Trees

This script reads one or more Newick phylogenetic trees, resolves multifurcations, and writes binary trees. This script is useful as a post-processing step when tree completion may introduce multifurcations. It resolves any internal node with more than two children and returns binary-refined Newick trees. Used for the evaluation purposes.

## Usage

```bash
python binarize_trees.py input_trees.nwk output_trees.nwk
````

## Input

The input file should contain one or more Newick trees, each ending with `;`.

Example:

```text
((A:1,B:1,C:1):1,D:2);
(A,B,C,D);
```

## Output

The output file contains the same trees with multifurcations resolved.

Example input:

```text
(A:1,B:1,C:1,D:1);
```

Example output:

```text
(D:1,((A:1,B:1):0,C:1):0);
```

The exact binary refinement is deterministic and based on leaf names.

## Notes

If the input tree has branch lengths, newly inserted internal branches get length `0`.

If the input tree has no branch lengths, newly inserted internal branches are written without branch lengths.

The script preserves existing leaf labels, internal labels, support values, and branch-length strings as plain Newick suffixes, but it does not preserve whitespace or line wrapping.

The added binary splits are artificial refinements of multifurcations. They should not be interpreted as biologically inferred relationships.
