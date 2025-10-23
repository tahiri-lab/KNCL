# Tree completion vs. tree pruning comparison

This folder contains scripts, datasets, and pre-computed pkl files (distance values and conflict tree pairs) for evaluating the agreement between BSD(*k*-NCL) (tree completion) and BSD(-) (tree pruning).

Three scenarios are evaluated based on how each distance interprets the similarity between two trees, $T_1$ and $T_2$, relative to a reference tree $T^*$ (the constructed supertree).

* **Scenario 1 (disagreement in ordering):**

  ```math
  \begin{aligned}
  \text{BSD}(k\text{-NCL})(T_1, T^*) &< \text{BSD}(k\text{-NCL})(T_2, T^*) \\
  \text{BSD(-)}(T_1, T^*) &> \text{BSD(-)}(T_2, T^*)
  \end{aligned}
  \quad \text{or} \quad
  \begin{aligned}
  \text{BSD}(k\text{-NCL})(T_2, T^*) &< \text{BSD}(k\text{-NCL})(T_1, T^*) \\
  \text{BSD(-)}(T_2, T^*) &> \text{BSD(-)}(T_1, T^*)
  \end{aligned}


* **Scenario 2 (different *k*-NCL distances, same BSD(-)):**

  ```math
  \begin{aligned}
  \text{BSD}(k\text{-NCL})(T_1, T^*) &\neq \text{BSD}(k\text{-NCL})(T_2, T^*) \\
  \text{BSD(-)}(T_1, T^*) &= \text{BSD(-)}(T_2, T^*)
  \end{aligned}

* **Scenario 3 (same *k*-NCL distance, different BSD(-)):**

  ```math
  \begin{aligned}
  \text{BSD}(k\text{-NCL})(T_1, T^*) &= \text{BSD}(k\text{-NCL})(T_2, T^*) \\
  \text{BSD(-)}(T_1, T^*) &\neq \text{BSD(-)}(T_2, T^*)
  \end{aligned}

The results include a table, line graphs, and violin charts showing the proportion of conflicts for each overlap level.
