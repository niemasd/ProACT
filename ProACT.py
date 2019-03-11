#!/usr/bin/env python3
from collections import deque
from copy import copy
from numpy import gradient
from treeswift import read_tree_newick,Node
from random import shuffle
from queue import PriorityQueue,Queue
from scipy.stats import norm
from sys import stderr
from warnings import warn
import heapq
INVALID_DATE = "Invalid date. Dates must be floats/integers"
NO_DIAG = "No diagnosis file specified, so using purely tree-based method"

# sort by edge length, then by difference with max infection time below sibling
def edgelength_diff_max_sib_inf(tree,inf,n):
    leaves = list()
    for u in tree.traverse_postorder():
        if u.is_leaf():
            u.max_inf = inf[u.label]; leaves.append(u)
        else:
            u.max_inf = max(c.max_inf for c in u.children)
    for u in leaves:
        u.max_sib_inf = max(c.max_inf for c in u.parent.children if c != u)
    score = {u.label: (u.edge_length, inf[u.label]-u.max_sib_inf) for u in leaves}
    return sorted(score.keys(), key=lambda x: score[x])[:n]

# sort by edge length, then by weighted/unweighted root-to-tip distance
def edgelength_norm_root_to_tip(tree,n):
    leaves = list()
    for u in tree.traverse_preorder():
        if u.is_root():
            u.root_dist = 0; u.root_dist_u = 0
        elif u.edge_length is None:
            u.root_dist = u.parent.root_dist; u.root_dist_u = u.parent.root_dist_u + 1
        else:
            u.root_dist = u.parent.root_dist + u.edge_length; u.root_dist_u = u.parent.root_dist_u + 1
        if u.is_leaf():
            leaves.append(u)
    score = {u.label: (u.edge_length, u.root_dist/u.root_dist_u) for u in leaves} # sort by edge length, then by weighted/unweighted root-to-tip distance
    return sorted(score.keys(), key=lambda x: score[x])[:n]

# run ProACT
if __name__ == "__main__":
    import argparse; from gzip import open as gopen
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, help="Input Tree File")
    parser.add_argument('-d', '--diagnosis', required=False, type=str, default=None, help="Input Diagnosis Time File")
    parser.add_argument('-n', '--number', required=False, type=str, default='All', help="Number of Individuals")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    args = parser.parse_args()
    if args.output == 'stdout':
        from sys import stdout; output = stdout
    else:
        output = open(args.output,'w')
    if args.diagnosis is None:
        inf = None; warn(NO_DIAG)
    else:
        if args.diagnosis.lower().endswith('.gz'):
            lines = [l.strip() for l in gopen(args.diagnosis).read().decode().strip().splitlines()]
        else:
            lines = [l.strip() for l in open(args.diagnosis).read().strip().splitlines()]
        inf = {}
        for line in lines:
            u,t = line.strip().split()
            try:
                inf[u] = float(t)
            except ValueError:
                raise ValueError(INVALID_DATE)
    tree = read_tree_newick(args.tree)
    num_leaves = len([l for l in tree.traverse_leaves()])
    if args.number == 'All':
        args.number = num_leaves
    else:
        args.number = int(args.number)
    assert 0 < args.number <= num_leaves, "Number of output individuals (%d) must be less than or equal to total number of individuals in tree (%d)" % (args.number,num_leaves)
    if inf is None:
        people = edgelength_norm_root_to_tip(tree,args.number)
    else:
        people = edgelength_diff_max_sib_inf(tree,inf,args.number)
    for u in people:
        output.write(u); output.write('\n')
