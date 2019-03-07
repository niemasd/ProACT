#!/usr/bin/env python3
'''
Find coefficients of linear combination of leaf properties for ProACT
'''
from common import individual_efficacy,individuals_from_tree,leaf_to_name,load_transmissions,optimal_order
from scipy.optimize import minimize
from sys import stderr
from treeswift import read_tree_newick
import argparse
VALS = ['edge_length','diff_root_to_tip']#,'diff_root_to_tip_u','num_sib_leaves']
TOP_N = 1000

def f(x, vals=None):
    score = {u:sum(x[i]*vals[VALS[i]][u] for i in range(len(x))) for u in vals['edge_length']}
    order = sorted(score.keys(), key=lambda u: score[u])
    tot_inf = sum(vals['infections'][u] for u in order[:TOP_N])
    if VERBOSE:
        stderr.write("Avg Infections: %f\n"%(tot_inf/TOP_N)); stderr.flush()
    return -tot_inf

def find_coefficients(tree,individuals,eff,maxit):
    assert maxit is None or maxit > 0, 'Max # of iterations must be positive'
    l2n = leaf_to_name(tree)
    vals = {v:dict() for v in ['infections','edge_length','diff_root_to_tip','diff_root_to_tip_u','num_sib_leaves']}
    for u in tree.traverse_postorder():
        if not u.is_root() and u.edge_length is None:
            u.edge_length = 0
        if u.is_leaf():
            u.num_leaves = 1
        else:
            u.num_leaves = sum(c.num_leaves for c in u.children)
    for u in tree.traverse_preorder():
        if u.is_root():
            u.root_dist = 0; u.root_dist_u = 0
        else:
            u.root_dist = u.parent.root_dist + u.edge_length; u.root_dist_u = u.parent.root_dist_u + 1
    avg_root_to_tip = sum(l.root_dist for l in tree.traverse_leaves())
    avg_root_to_tip_u = sum(l.root_dist_u for l in tree.traverse_leaves())
    for u in tree.traverse_leaves():
        uname = l2n[u]
        vals['edge_length'][uname] = u.edge_length
        vals['diff_root_to_tip'][uname] = abs(u.root_dist - avg_root_to_tip)
        vals['diff_root_to_tip_u'][uname] = abs(u.root_dist_u - avg_root_to_tip_u)
        vals['num_sib_leaves'][uname] = sum(c.num_leaves for c in u.parent.children if c != u)
        vals['infections'][uname] = eff[uname]
    x0 = [1.]*len(VALS)
    if maxit is None:
        result = minimize(f, x0, args=vals, method='Powell')
    else:
        result = minimize(f, x0, args=vals, method='Powell', options={'maxiter':maxit})
    return {VALS[i]:result.x[i] for i in range(len(VALS))}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tn', '--transmissions', required=True, type=str, help="Transmission Network (FAVITES format)")
    parser.add_argument('-tr', '--tree', required=True, type=str, help="Phylogenetic Tree (Newick Format)")
    parser.add_argument('-t', '--from_time', required=True, type=float, help="From Time")
    parser.add_argument('-tt', '--to_time', required=False, type=float, default=float('inf'), help="To Time")
    parser.add_argument('-m', '--max_it', required=False, type=int, default=None, help="Max # of Iterations")
    parser.add_argument('-v', '--verbose', action='store_true', help="Verbose Mode")
    args = parser.parse_args()
    global VERBOSE; VERBOSE = args.verbose
    tree = read_tree_newick(args.tree)
    individuals = individuals_from_tree(tree)
    trans = load_transmissions(args.transmissions)
    eff = individual_efficacy(individuals,trans,args.from_time,args.to_time)
    for pair in find_coefficients(tree,individuals,eff,args.max_it).items():
        print('%s\t%f'%pair)
