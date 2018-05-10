#!/usr/bin/env python3
'''
Given a tree T, a start time S, and an end time E, remove all leaves from T that
are outside of the time window [S,E]. Trees must be in the Newick format, and
leaf labels must be delimited by the '|' character, with the last field being
the time associated with the leaf.
'''
from dendropy import Node,Tree
from queue import Queue

# get the label of Dendropy node u
def label(u):
    return str(u.taxon).replace("'",'')

# my own version of extract_tree_with_taxa
def extract_tree_with_taxa(tree, taxa, suppress_unifurcations=True):
    taxon_to_leaf = {}
    for n in tree.preorder_node_iter():
        n.keep = False
        if n.is_leaf():
            taxon_to_leaf[n.taxon] = n
    for t in taxa:
        for n in taxon_to_leaf[t].ancestor_iter(inclusive=True):
            n.keep = True
    out = Tree()
    q_old = Queue(); q_old.put(tree.seed_node)
    q_new = Queue(); q_new.put(out.seed_node)
    while not q_old.empty():
        n_old = q_old.get(); n_new = q_new.get()
        for c_old in n_old.child_node_iter():
            if c_old.keep:
                c_new = Node(taxon=c_old.taxon, label=c_old.label, edge_length=c_old.edge_length); n_new.add_child(c_new)
                q_old.put(c_old); q_new.put(c_new)
    if suppress_unifurcations:
        out.suppress_unifurcations()
    return out

# run main program
if __name__ == "__main__":
    import argparse; from gzip import open as gopen
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str, default='stdin', help="Input Tree File")
    parser.add_argument('-s', '--start', required=False, type=float, default=float('-inf'), help="Window Start Time")
    parser.add_argument('-e', '--end', required=False, type=float, default=float('inf'), help="Window End Time")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    args = parser.parse_args()
    assert args.start != float('-inf') or args.end != float('inf'), "Must specify either start time or end time (or both)"
    if args.tree.lower() == 'stdin':
        from sys import stdin; tree = Tree.get(data=stdin.read(),schema='newick')
    elif args.tree.lower().endswith('.gz'):
        tree = Tree.get(data=gopen(args.tree).read().decode(),schema='newick')
    else:
        tree = Tree.get(data=open(args.tree).read(),schema='newick')
    if args.output == 'stdout':
        from sys import stdout; output = stdout
    else:
        output = open(args.output,'w')
    taxa = [leaf.taxon for leaf in tree.leaf_node_iter() if args.start < float(label(leaf).split('|')[-1]) < args.end]
    #output.write((tree.extract_tree_with_taxa(taxa, suppress_unifurcations=False).as_string('newick')))
    output.write(extract_tree_with_taxa(tree,taxa).as_string('newick').replace("'",''))
    output.write('\n'); output.close()