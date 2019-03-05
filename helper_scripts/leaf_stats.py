#!/usr/bin/env python3
from common import individual_efficacy,leaf_to_name,load_transmissions
from treeswift import read_tree_newick
import argparse
ORDER = ['name', 'efficacy', 'edge_length', 'root_to_tip', 'root_to_tip_u', 'sib_leaves']

def compute_sib_leaves(tree):
    num_leaves = dict()
    for u in tree.traverse_postorder():
        if u.is_leaf():
            num_leaves[u] = 1
        else:
            num_leaves[u] = sum(num_leaves[c] for c in u.children)
    return {L2N[l]:sum(num_leaves[c] for c in l.parent.children if c != l) for l in L2N}

def compute_root_to_tip(tree, weighted=True):
    root_to_tip = dict()
    for u in tree.traverse_preorder():
        if u.is_root():
            root_to_tip[u] = 0
        elif not weighted:
            root_to_tip[u] = root_to_tip[u.parent] + 1
        elif u.edge_length is None:
            root_to_tip[u] = root_to_tip[u.parent]
        else:
            root_to_tip[u] = root_to_tip[u.parent] + u.edge_length
    return {L2N[l]:root_to_tip[l] for l in L2N}

def compute_vals(tree):
    vals = dict()
    vals['sib_leaves'] = compute_sib_leaves(tree)
    vals['edge_length'] = {L2N[l]:l.edge_length for l in L2N}
    vals['root_to_tip'] = compute_root_to_tip(tree, weighted=True)
    vals['root_to_tip_u'] = compute_root_to_tip(tree, weighted=False)
    return vals

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tr', '--tree', required=True, type=str, help="Phylogenetic Tree (Newick format)")
    parser.add_argument('-tn', '--transmissions', required=True, type=str, help="Transmission Network (FAVITES format)")
    parser.add_argument('-t', '--from_time', required=True, type=float, help="From Time (for # transmissions)")
    parser.add_argument('-tt', '--to_time', required=False, type=float, default=float('inf'), help="To Time (for # transmissions)")
    args = parser.parse_args()
    tree = read_tree_newick(args.tree)
    global L2N; L2N = leaf_to_name(tree)
    vals = compute_vals(tree)
    vals['name'] = {L2N[l]:L2N[l] for l in L2N}
    vals['efficacy'] = individual_efficacy([L2N[l] for l in tree.traverse_leaves()],load_transmissions(args.transmissions),args.from_time,args.to_time)
    print(','.join(ORDER))
    for l in L2N:
        ln = L2N[l]; print('%s,%d,%f,%f,%d,%d' % tuple(vals[k][ln] for k in ORDER))
