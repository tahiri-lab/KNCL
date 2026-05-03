#!/usr/bin/env python3

"""
Usage:
    python binarize_trees.py input_trees.nwk output_trees.nwk

Reads one or more Newick trees from input_trees.nwk, resolves any
multifurcations, and writes binary-refined Newick trees to output_trees.nwk.

New internal branches are given length 0 if the input tree has branch lengths.
Otherwise, no branch length is added to new internal branches.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field

@dataclass
class Node:
    children: list["Node"] = field(default_factory=list)
    suffix: str = ""

    def is_leaf(self) -> bool:
        return len(self.children) == 0

def split_newick_trees(text: str) -> list[str]:
    """
    Split text into Newick trees using semicolons outside quotes/comments.
    """
    trees = []
    buf = []

    in_quote = False
    in_comment = False
    i = 0

    while i < len(text):
        ch = text[i]
        buf.append(ch)

        if in_comment:
            if ch == "]":
                in_comment = False

        elif in_quote:
            if ch == "'":
                # Escaped quote inside a quoted label: ''
                if i + 1 < len(text) and text[i + 1] == "'":
                    i += 1
                    buf.append(text[i])
                else:
                    in_quote = False

        else:
            if ch == "[":
                in_comment = True
            elif ch == "'":
                in_quote = True
            elif ch == ";":
                tree = "".join(buf).strip()
                if tree:
                    trees.append(tree)
                buf = []

        i += 1

    leftover = "".join(buf).strip()
    if leftover:
        raise ValueError("Input contains a tree without a terminating semicolon.")

    return trees

def has_branch_lengths(newick: str) -> bool:
    """
    Return True if the Newick string contains ':' outside quotes/comments.
    """
    in_quote = False
    in_comment = False
    i = 0

    while i < len(newick):
        ch = newick[i]

        if in_comment:
            if ch == "]":
                in_comment = False

        elif in_quote:
            if ch == "'":
                if i + 1 < len(newick) and newick[i + 1] == "'":
                    i += 1
                else:
                    in_quote = False

        else:
            if ch == "[":
                in_comment = True
            elif ch == "'":
                in_quote = True
            elif ch == ":":
                return True

        i += 1

    return False

class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.pos = 0

    def parse(self) -> Node:
        node = self.parse_subtree()
        self.skip_spaces()

        if self.pos >= len(self.text) or self.text[self.pos] != ";":
            raise ValueError("Expected semicolon at end of tree.")

        return node

    def parse_subtree(self) -> Node:
        self.skip_spaces()

        if self.pos >= len(self.text):
            raise ValueError("Unexpected end of Newick string.")

        if self.text[self.pos] == "(":
            self.pos += 1
            children = []

            while True:
                child = self.parse_subtree()
                children.append(child)

                self.skip_spaces()

                if self.pos >= len(self.text):
                    raise ValueError("Unexpected end inside internal node.")

                if self.text[self.pos] == ",":
                    self.pos += 1
                    continue

                if self.text[self.pos] == ")":
                    self.pos += 1
                    break

                raise ValueError(f"Unexpected character: {self.text[self.pos]!r}")

            suffix = self.parse_suffix()
            return Node(children=children, suffix=suffix)

        else:
            suffix = self.parse_suffix()
            if not suffix:
                raise ValueError("Found an empty leaf label.")
            return Node(children=[], suffix=suffix)

    def parse_suffix(self) -> str:
        """
        Parse node label / support / branch length until ',', ')' or ';'.
        """
        start = self.pos
        in_quote = False
        in_comment = False

        while self.pos < len(self.text):
            ch = self.text[self.pos]

            if in_comment:
                if ch == "]":
                    in_comment = False

            elif in_quote:
                if ch == "'":
                    if self.pos + 1 < len(self.text) and self.text[self.pos + 1] == "'":
                        self.pos += 1
                    else:
                        in_quote = False

            else:
                if ch == "[":
                    in_comment = True
                elif ch == "'":
                    in_quote = True
                elif ch in ",);":
                    break

            self.pos += 1

        return self.text[start:self.pos].strip()

    def skip_spaces(self):
        while self.pos < len(self.text) and self.text[self.pos].isspace():
            self.pos += 1


def leaf_name_from_suffix(suffix: str) -> str:
    """
    Extract the leaf name from a leaf suffix such as:
        A:0.1
        'A species':0.1
        A[comment]:0.1
    """
    suffix = suffix.strip()

    if not suffix:
        raise ValueError("Unnamed leaf found.")

    if suffix.startswith("'"):
        name_chars = []
        i = 1

        while i < len(suffix):
            ch = suffix[i]

            if ch == "'":
                if i + 1 < len(suffix) and suffix[i + 1] == "'":
                    name_chars.append("'")
                    i += 2
                    continue
                return "".join(name_chars)

            name_chars.append(ch)
            i += 1

        raise ValueError("Unterminated quoted leaf name.")

    name_chars = []

    for ch in suffix:
        if ch in ":[":
            break
        name_chars.append(ch)

    name = "".join(name_chars).strip()

    if not name:
        raise ValueError("Unnamed leaf found.")

    return name


def subtree_leaf_key(node: Node) -> tuple[str, ...]:
    """
    Deterministic key based on the sorted leaf names below a node.
    """
    if node.is_leaf():
        return (leaf_name_from_suffix(node.suffix),)

    names = []
    for child in node.children:
        names.extend(subtree_leaf_key(child))

    return tuple(sorted(names))


def binarize(node: Node, new_internal_suffix: str) -> int:
    """
    Resolve all multifurcations below this node.

    Returns the number of newly inserted internal nodes.
    """
    inserted = 0

    for child in node.children:
        inserted += binarize(child, new_internal_suffix)

    while len(node.children) > 2:
        children_sorted = sorted(node.children, key=subtree_leaf_key)

        c1 = children_sorted[0]
        c2 = children_sorted[1]

        remaining = [c for c in node.children if c is not c1 and c is not c2]

        new_node = Node(
            children=[c1, c2],
            suffix=new_internal_suffix,
        )

        node.children = remaining + [new_node]
        inserted += 1

    return inserted


def write_newick(node: Node) -> str:
    if node.is_leaf():
        return node.suffix

    return "(" + ",".join(write_newick(child) for child in node.children) + ")" + node.suffix


def count_multifurcations(node: Node) -> int:
    count = 1 if len(node.children) > 2 else 0

    for child in node.children:
        count += count_multifurcations(child)

    return count


def main() -> int:
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage:\n"
            "    python binarize_trees.py input_trees.nwk output_trees.nwk\n"
        )
        return 1

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(input_file, "r", encoding="utf-8") as f:
        text = f.read()

    try:
        tree_strings = split_newick_trees(text)
    except Exception as exc:
        sys.stderr.write(f"ERROR while reading input trees: {exc}\n")
        return 1

    if not tree_strings:
        sys.stderr.write("ERROR: no Newick trees found.\n")
        return 1

    output_trees = []

    for i, tree_string in enumerate(tree_strings, start=1):
        try:
            parser = NewickParser(tree_string)
            tree = parser.parse()

            before = count_multifurcations(tree)

            # If the input tree uses branch lengths, newly inserted internal
            # branches get length 0. Otherwise, they are written without length.
            new_internal_suffix = ":0" if has_branch_lengths(tree_string) else ""

            inserted = binarize(tree, new_internal_suffix)

            after = count_multifurcations(tree)

            if after != 0:
                raise RuntimeError("Tree still contains multifurcations after binarization.")

            output_trees.append(write_newick(tree) + ";")

            sys.stderr.write(
                f"Tree {i}: multifurcations={before}, "
                f"inserted_internal_nodes={inserted}\n"
            )

        except Exception as exc:
            sys.stderr.write(f"ERROR while processing tree {i}: {exc}\n")
            return 1

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("\n".join(output_trees) + "\n")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
