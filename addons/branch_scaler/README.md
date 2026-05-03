# Branch-Length Scaling for *k*-NCL Tree Completion

`branch_scaler.py` is an optional preprocessing and postprocessing utility for *k*-NCL-based tree completion.

It is useful when two input phylogenetic trees have branch lengths on very different numerical scales. The script can scale the two original input trees before completion, save the scale and unscale factors in a JSON file, and later use that JSON file to return the completed trees to their original branch-length scales.
