#!/usr/bin/env python3

import os, sys
from ete3 import Tree


def main():
    nwk_file = sys.argv[1]
    tissues_file = sys.argv[2]
    out_file = sys.argv[3]

    tree = Tree(nwk_file, format=3)

    tissues = {}
    with open(tissues_file, "r") as file:
        for line in file:
            fields = line.strip().split(" ")
            tissues[fields[0]] = fields[1]

    for node in tree.traverse():
        name = node.name
        node.name = name + "_" + tissues[name]

    tree.write(outfile=out_file, format=8)


if __name__ == "__main__":
    main()
