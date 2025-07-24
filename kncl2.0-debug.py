# Debug version
"""
k-Nearest Common Leaves (k-NCL) Tree Completion Algorithm

This script implements the k-NCL algorithm for completing phylogenetic trees
that are defined on different but overlapping taxon sets. The implementation
uses the ete3 library to work with phylogenetic trees.
"""

from functools import lru_cache
from math import floor, log2
from ete3 import Tree

# ──────────────────────────────────────────────────────────────────────
#  Helper functions
# ──────────────────────────────────────────────────────────────────────
def _refresh_name_cache(tree: Tree):
    tree._nmap = {n.name: n for n in tree.traverse() if n.name}


def ensure_unique_internal_node_names(tree: Tree, prefix="IN"):
    counter = 1
    for node in tree.traverse("postorder"):
        if not node.is_leaf() and (not node.name or node.name.startswith(prefix)):
            node.name = f"{prefix}_{counter}"
            counter += 1
    _refresh_name_cache(tree)


def find_common_leaves(T1: Tree, T2: Tree):
    return set(T1.get_leaf_names()).intersection(T2.get_leaf_names())


def find_distinct_leaves(T: Tree, common):
    return set(T.get_leaf_names()) - common


def _is_descendant(node, ancestor):
    return node == ancestor or ancestor in node.iter_ancestors()


# ──────────────────────────────────────────────────────────────────────
#  Find maximal distinct‑leaf subtrees
# ──────────────────────────────────────────────────────────────────────
def findSD(tree: Tree, dl_set):
    for n in tree.traverse("postorder"):
        if n.is_leaf():
            n.allDistinct = n.name in dl_set
        else:
            n.allDistinct = all(ch.allDistinct for ch in n.children)
    return [n for n in tree.traverse("postorder")
        if n.allDistinct and (n.up is None or not n.up.allDistinct)]


# ──────────────────────────────────────────────────────────────────────
#  Depth & LCA oracle
# ──────────────────────────────────────────────────────────────────────
class DistOracle:
    __slots__ = ("dist_to_root", "_euler", "_depth", "_first", "_st", "_log2", "_dist_cached")

    def __init__(self, tree: Tree, cache_size: int = 20_000):
        self.dist_to_root = {}
        self._annotate_depth(tree)
        self._build_lca_struct(tree)

        @lru_cache(maxsize=cache_size)
        def _dist_cached(a, b):
            lca = self._lca(a, b)
            return self.dist_to_root[a] + self.dist_to_root[b] - 2 * self.dist_to_root[lca]

        self._dist_cached = _dist_cached

    def dist(self, a, b):
        return self._dist_cached(a, b)

    def _annotate_depth(self, tree):
        for n in tree.traverse("preorder"):
            self.dist_to_root[n] = 0.0 if n.up is None else self.dist_to_root[n.up] + n.dist

    def _build_lca_struct(self, tree):
        euler, depth, first = [], [], {}

        def dfs(v, d):
            first.setdefault(v, len(euler))
            euler.append(v)
            depth.append(d)
            for ch in v.children:
                dfs(ch, d + 1)
                euler.append(v)
                depth.append(d)

        dfs(tree, 0)
        self._euler, self._depth, self._first = euler, depth, first

        m = len(euler)
        kmax = floor(log2(m)) + 1
        st = [[0] * m for _ in range(kmax)]
        st[0] = list(range(m))
        for k in range(1, kmax):
            half = 1 << (k - 1)
            for i in range(m - (1 << k) + 1):
                a, b = st[k - 1][i], st[k - 1][i + half]
                st[k][i] = a if depth[a] < depth[b] else b
        self._st = st
        self._log2 = [0] * (m + 1)
        for i in range(2, m + 1):
            self._log2[i] = self._log2[i >> 1] + 1

    def _lca(self, a, b):
        ia, ib = self._first[a], self._first[b]
        if ia > ib:
            ia, ib = ib, ia
        span = ib - ia + 1
        k = self._log2[span]
        left = self._st[k][ia]
        right = self._st[k][ib - (1 << k) + 1]
        return self._euler[left] if self._depth[left] < self._depth[right] else self._euler[right]

    def dist_leaf_to_node(self, leaf, node):
        return self.dist(leaf, node)


# ──────────────────────────────────────────────────────────────────────
#  Adjustment rate functions
# ──────────────────────────────────────────────────────────────────────
def compute_global_adjustment_rate(T1, T2, common_names, D1, D2):
    L1 = [T1 & n for n in common_names]
    L2 = [T2 & n for n in common_names]
    num = den = 0.0
    for i in range(len(L1)):
        for j in range(i + 1, len(L1)):
            num += D1.dist(L1[i], L1[j])
            den += D2.dist(L2[i], L2[j])
    return 1.0 if abs(den) < 1e-15 else num / den


def compute_leaf_based_adjustment_rate(Dtgt, Dsrc, leaf_tgt, leaf_src, common_names, tree_tgt, tree_src):
    num = den = 0.0
    for name in common_names:
        if name == leaf_tgt.name:
            continue
        num += Dtgt.dist(leaf_tgt, tree_tgt & name)
        den += Dsrc.dist(leaf_src, tree_src & name)
    return 1.0 if abs(den) < 1e-15 else num / den


# ──────────────────────────────────────────────────────────────────────
#  Misc helpers
# ──────────────────────────────────────────────────────────────────────
def scale_subtree(root, factor):
    for n in root.traverse():
        n.dist *= factor


def get_k_nearest_common_leaves(Dsrc, root, common_names, tree_src, k):
    pairs = [(name, Dsrc.dist_leaf_to_node(tree_src & name, root)) for name in common_names]
    pairs.sort(key=lambda t: t[1])
    if not pairs:
        return []
    kth = pairs[min(k, len(pairs)) - 1][1]
    return [n for n, d in pairs if d <= kth + 1e-15]


# ──────────────────────────────────────────────────────────────────────
#  Insertion & Optimization
# ──────────────────────────────────────────────────────────────────────
def find_optimal_insertion_point(target_tree, Dtgt, subtree_root, ncl_names, d_p, dl_names, min_terminal=1e-3):
    best_edge, best_x, best_val = None, 0.0, float("inf")
    dl_names = set(dl_names)
    m = len(ncl_names)
    if m == 0:
        return None, 0.0, float("inf")

    for child in target_tree.traverse("postorder"):
        parent = child.up
        if parent is None:
            continue
        if getattr(child, "inInserted", False):
            continue
        if set(child.get_leaf_names()) & dl_names:
            continue

        e_len = child.dist
        if e_len < 1e-15:
            continue

        numerator = 0.0
        for name in ncl_names:
            leaf = target_tree & name
            A_c = Dtgt.dist_leaf_to_node(leaf, parent)
            eps = -1 if _is_descendant(leaf, child) else 1
            numerator += eps * (d_p[name] - A_c)
        x_opt = numerator / (e_len * m)

        if child.is_leaf():
            x_opt = max(0.0, min(x_opt, 1 - min_terminal / e_len))
        else:
            x_opt = max(0.0, min(x_opt, 1.0))

        for x in (x_opt, 0.0):
            of_val = 0.0
            for name in ncl_names:
                leaf = target_tree & name
                A_c = Dtgt.dist_leaf_to_node(leaf, parent)
                eps = -1 if _is_descendant(leaf, child) else 1
                obs = A_c + eps * x * e_len
                diff = obs - d_p[name]
                of_val += diff * diff
            if of_val < best_val - 1e-12:
                best_val, best_edge, best_x = of_val, (parent, child), x

    return best_edge, best_x, best_val


def insert_subtree_at_point(target_tree, edge, x_opt, subtree_copy):
    parent, child = edge
    orig_len = child.dist
    d_up, d_down = x_opt * orig_len, (1 - x_opt) * orig_len

    if x_opt < 1e-15:
        parent.add_child(subtree_copy)
    else:
        mid = Tree()
        mid.dist = d_up
        parent.remove_child(child)
        parent.add_child(mid)
        child.dist = d_down
        mid.add_child(child)
        mid.add_child(subtree_copy)

    _refresh_name_cache(target_tree)


# ──────────────────────────────────────────────────────────────────────
#  Main k-NCL Function
# ──────────────────────────────────────────────────────────────────────
def kNCL(T1_raw: Tree, T2_raw: Tree, k=None):
    T1, T2 = T1_raw.copy(), T2_raw.copy()
    ensure_unique_internal_node_names(T1, "T1IN")
    ensure_unique_internal_node_names(T2, "T2IN")

    common = find_common_leaves(T1, T2)
    print("Common leaves:", sorted(common))
    if len(common) < 2:
        raise ValueError("Need at least two common leaves")

    if k is None:
        k = len(common)
    if k < 2:
        raise ValueError("k must be ≥ 2")

    DL1, DL2 = find_distinct_leaves(T1, common), find_distinct_leaves(T2, common)
    SD1, SD2 = findSD(T1, DL1), findSD(T2, DL2)

    print("\nMaximal distinct-leaf subtrees in T1:")
    for sd in SD1:
        print("  -", sd.get_leaf_names())

    print("\nMaximal distinct-leaf subtrees in T2:")
    for sd in SD2:
        print("  -", sd.get_leaf_names())

    D1, D2 = DistOracle(T1), DistOracle(T2)
    r12 = compute_global_adjustment_rate(T1, T2, common, D1, D2)
    r21 = compute_global_adjustment_rate(T2, T1, common, D2, D1)

    print(f"\nGlobal adjustment rate T1→T2: {r12:.4f}")
    print(f"Global adjustment rate T2→T1: {r21:.4f}")

    def _insert_all(target_tree, source_tree, subtrees, global_r, Dtgt, Dsrc, direction):
        print(f"\n[{direction}] Global adjustment rate = {global_r:.4f}")
        for idx, S in enumerate(subtrees):
            print(f"\n[{direction}] Processing subtree {idx + 1}/{len(subtrees)} with leaves: {S.get_leaf_names()}")
            S_copy = S.copy()
            scale_subtree(S_copy, global_r)
            for n in S_copy.traverse():
                n.add_feature("inInserted", True)

            ncl = get_k_nearest_common_leaves(Dsrc, S, common, source_tree, k)
            print(f"[{direction}] Nearest common leaves: {ncl}")
            if not ncl:
                continue

            att_node = S.up if S.up else S
            d_p = {}
            for name in ncl:
                leaf_src = source_tree & name
                leaf_tgt = target_tree & name
                r_lc = compute_leaf_based_adjustment_rate(
                    Dtgt, Dsrc, leaf_tgt, leaf_src, common, target_tree, source_tree
                )
                d = Dsrc.dist_leaf_to_node(leaf_src, att_node) * r_lc
                d_p[name] = d
                print(f"[{direction}] Leaf {name}: rate = {r_lc:.4f}, adjusted dist = {d:.4f}")

            edge, x_opt, of_val = find_optimal_insertion_point(
                target_tree, Dtgt, S, ncl, d_p, find_distinct_leaves(target_tree, common)
            )
            if edge:
                print(f"[{direction}] Inserting at edge ({edge[0].name}, {edge[1].name}) with x = {x_opt:.4f}")
                insert_subtree_at_point(target_tree, edge, x_opt, S_copy)
                Dtgt.__init__(target_tree)
                _refresh_name_cache(target_tree)
            else:
                print(f"[{direction}] No valid insertion point found.")

    _insert_all(T1, T2, SD2, r12, D1, D2, "T1 ← T2")
    _insert_all(T2, T1, SD1, r21, D2, D1, "T2 ← T1")

    for t in (T1, T2):
        for n in t.traverse():
            if not n.is_leaf():
                n.name = ""

    return T1, T2

def run_kNCL_on_newick_strings(newick1, newick2, k):
    T1 = Tree(newick1, format=1)
    T2 = Tree(newick2, format=1)
    T1_plus, T2_plus = kNCL(T1, T2, k)
    print("\nCompleted Tree 1:")
    print(T1_plus.write(format=1))
    print("\nCompleted Tree 2:")
    print(T2_plus.write(format=1))
