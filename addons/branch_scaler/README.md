# Branch-length scaling for *k*-NCL tree completion

`branch_scaler.py` is an optional preprocessing and postprocessing utility for *k*-NCL-style tree completion.

It is useful when a collection of input phylogenetic trees has branch lengths on very different numerical scales. The script scales each tree independently, writes the scaled trees to a new Newick file, and saves the corresponding scale/unscale factors in a JSON file.

The scaled trees can then be used directly as input to *k*-NCL. After tree completion, the completed trees can be returned to their original branch-length scales using the saved JSON factor file.

## Basic idea

The script applies multiplicative branch-length scaling.

For each tree:

```text
scaled_branch_length = original_branch_length * scale_factor
````

Later, to return to the original scale:

```text
original_scale_branch_length = scaled_branch_length * unscale_factor
```

where:

```text
unscale_factor = 1 / scale_factor
```

This script does not change topology. It only multiplies branch lengths.

## Requirements

Python 3. No external Python packages are required.

The input file should contain one or more Newick trees. Each tree must end with a semicolon.

Example:

```text
((A:1,B:1):1,C:2);
((A:10,B:10):10,C:20);
(A:1,B:2,C:3,D:4);
```

All non-root branches must have branch lengths. Negative branch lengths are not allowed.

## Scaling methods

The script supports three independent scaling methods.

### 1. `median-all-pairwise`

This is the recommended default.

For each tree, the script computes all pairwise leaf-to-leaf distances and scales the tree by:

```text
scale_factor = 1 / median(all pairwise leaf-to-leaf distances)
```

<!-- After scaling, the median pairwise leaf distance of each tree becomes 1. -->

Use this method when you want a robust default. It is less sensitive to outlier taxa and unusually long branches than the mean.

Recommended for most *k*-NCL preprocessing workflows.

### 2. `mean-all-pairwise`

For each tree, the script computes all pairwise leaf-to-leaf distances and scales the tree by:

```text
scale_factor = 1 / mean(all pairwise leaf-to-leaf distances)
```

<!-- After scaling, the mean pairwise leaf distance of each tree becomes 1. -->

Use this method when branch lengths are relatively homogeneous and you do not expect strong outliers.

This method uses the whole distance structure of the tree, but it is more sensitive to long branches than `median-all-pairwise`.

### 3. `total-tree-length`

For each tree, the script computes the total non-root branch length and scales the tree by:

```text
scale_factor = 1 / total_tree_length
```

<!-- After scaling, the total tree length becomes 1. -->

This method is simple and fast. However, it is more sensitive to the number of taxa and the amount of tree resolution. Trees with more taxa usually have more branches, so total tree length may reflect taxon sampling as well as branch-length scale.

Use this method when specifically total-length normalization is needed.

## Recommendations

Use `median-all-pairwise` as the default choice.

Use `mean-all-pairwise` when the trees have no extreme branch-length outliers and you want the average leaf-to-leaf distance to be normalized.

Use `total-tree-length` when you want a simple whole-tree normalization, but be cautious if trees differ strongly in taxon number or sampling density.

See more detailed recommendations on provided scaling methods below (at the end).

Do not use branch-length scaling to force together trees whose branch lengths have incompatible meanings. For example, if one tree has branch lengths in substitutions/site and another tree has branch lengths in absolute time, a simple multiplicative scaling may not be biologically meaningful.

## Step 1: Scale all trees in a file

Example:

```bash
python branch_scaler.py scale \
  --input input_trees.nwk \
  --output input_trees.scaled.nwk \
  --factors scaling_factors.json \
  --method median-all-pairwise
```

This writes:

```text
input_trees.scaled.nwk
scaling_factors.json
```

The scaled Newick file contains the same number of trees as the input file, in the same order.

The factor file records one scale factor and one unscale factor per tree.

## Example factor file

The JSON factor file has entries like this:

```json
{
  "version": 2,
  "mode": "independent-batch",
  "method": "median-all-pairwise",
  "meaning": "scaled_branch_length = original_branch_length * scale_factor",
  "tree_count": 3,
  "trees": [
    {
      "tree_index": 1,
      "scale_factor": 0.25,
      "unscale_factor": 4.0
    },
    {
      "tree_index": 2,
      "scale_factor": 0.025,
      "unscale_factor": 40.0
    },
    {
      "tree_index": 3,
      "scale_factor": 0.2,
      "unscale_factor": 5.0
    }
  ]
}
```

Tree indices are 1-based.

## Step 2: Run *k*-NCL on scaled trees

Use the scaled trees as *k*-NCL input.

Suppose *k*-NCL produces:

```text
completed_scaled.nwk
```

containing the completed versions of scaled tree, in the original order.

## Step 3: Unscale completed trees

Use the saved factor file.

If the completed file contains the same number of trees in the same order as the original scaled file, you can run:

```bash
python branch_scaler.py unscale \
  --input completed_scaled.nwk \
  --output completed_unscaled.nwk \
  --factors scaling_factors.json
```

If the completed file contains only selected trees, provide their original tree indices:

```bash
python branch_scaler.py unscale \
  --input completed_scaled.nwk \
  --output completed_unscaled.nwk \
  --factors scaling_factors.json \
  --tree-indices 3,7
```

This means:

```text
first tree in completed_scaled.nwk  -> use the unscale factor from original tree 3
second tree in completed_scaled.nwk -> use the unscale factor from original tree 7
```

The output file will contain the completed trees returned to their original respective branch-length scales.

## Full recommended workflow

```text
1. Start with a file containing multiple original Newick trees.

2. Scale all trees independently:

   python branch_scaler.py scale \
     --input input_trees.nwk \
     --output input_trees.scaled.nwk \
     --factors scaling_factors.json \
     --method median-all-pairwise

3. Use selected scaled trees as input to k-NCL.

4. Unscale the completed k-NCL outputs using the saved factors:

   python branch_scaler.py unscale \
     --input completed_scaled.nwk \
     --output completed_unscaled.nwk \
     --factors scaling_factors.json \
     --tree-indices 3,7
```

## Additional notes

* The script preserves tree topology and leaf labels.

* The script rewrites Newick output, so exact whitespace and line wrapping are not preserved. It means the script preserves the tree content, but not the exact text formatting of the original file.

* All non-root branches must have branch lengths during scaling.

* The script scales every branch length it finds during both scaling and unscaling. Therefore, completed *k*-NCL trees should include branch lengths on inserted branches too.

* If a tree has a zero median pairwise distance, `median-all-pairwise` cannot be used for that tree. In that case, try `mean-all-pairwise` or `total-tree-length` if biologically appropriate.

* The added or completed branch lengths after *k*-NCL should be interpreted with care. Scaling makes branch-length magnitudes more comparable numerically, but it does not guarantee that branch lengths from different biological sources have the same meaning.

# Empirical recommendations on branch length scaling

## When is extra scaling needed before *k*-NCL?

*k*-NCL already includes internal branch-length adjustment. In particular, it estimates scale relationships from distances among common leaves and uses those rates during subtree insertion. Therefore, external preprocessing is not always necessary.

External scaling is most useful when the input trees have large global branch-length differences or when the overlap between a pair of trees is small or potentially unrepresentative.

### Cases where internal *k*-NCL scaling may be enough

Extra scaling is usually optional when two trees have:

* the same branch-length meaning
* enough common taxa (see details below)
* similar relative distance structure among common taxa
* a mostly multiplicative scale difference

For example, if most common-leaf distances satisfy approximately:

```text
distance_tree_A(i,j) ≈ c * distance_tree_B(i,j)
```

for a relatively stable constant `c`, then *k*-NCL's internal adjustment may already handle the scale difference reasonably well.

#### What counts as enough common taxa?

As a rough guideline:

- 2 common taxa gives only 1 pairwise distance and is not enough for reliable scale estimation.
- 3-5 common taxa is very weak and can be strongly affected by one unusual branch.
- 6-10 common taxa may be usable if the taxa are spread across the tree.
- 10-20 common taxa is usually more reliable.
- More than 20 common taxa is generally good, provided the shared taxa are not all concentrated in one small clade.

The shared taxa should ideally be distributed across the tree. A small number of well-distributed common taxa can be more useful than many common taxa restricted to one shallow clade.


### Cases where external scaling is recommended

External scaling is recommended when:

* tree branch lengths differ by more than about 10x
* the two trees share only a small number of leaves
* the common leaves are concentrated in one clade and may not represent the full tree
* branch lengths have very different numerical magnitudes, such as 0.0001 in one tree and 100 in another

### Practical diagnostic before scaling

Before deciding whether scaling is needed, compare simple branch-length summaries across trees:

* median all-leaf pairwise distance
* mean all-leaf pairwise distance
* total tree length

For two trees, compute ratios such as:

* median_distance_tree_A / median_distance_tree_B
* mean_distance_tree_A / mean_distance_tree_B
* total_length_tree_A / total_length_tree_B

As a practical rule of thumb (empirical):

```text
scale ratio < 2-5x:
    external scaling is usually optional

scale ratio around 5-10x:
    external scaling is reasonable, especially in batch tree processing

scale ratio > 10x:
    external scaling is recommended

scale ratio varies strongly depending on which summary is used:
    scaling may still help numerically, but completed branch lengths should be interpreted cautiously
```

These thresholds are practical guidelines, not strict mathematical rules.

### Check common-leaf distance compatibility

For a specific pair of trees, it can also be useful to compare distances among shared leaves.

For each pair of common leaves `(i,j)`, compute:

* distance_tree_A(i,j)
* distance_tree_B(i,j)
* ratio_ij = distance_tree_A(i,j) / distance_tree_B(i,j)

If the ratios `ratio_ij` are fairly consistent, then the scale mismatch is mostly global and scaling is meaningful.

If the ratios vary strongly, then the difference between the trees is not only a global scale issue. External scaling may still improve numerical comparability, but it cannot fully correct local rate heterogeneity, different branch-length estimation behavior, or conflicting biological signals.

### Important note

External scaling should not be used to force together branch lengths with incompatible meanings. For example:

* Tree 1 branch lengths = substitutions/site
* Tree 2 branch lengths = absolute time

A single multiplicative scale factor may not make such trees biologically comparable.

In such cases, the completed topology may still be useful, but completed branch lengths should be interpreted with care.
