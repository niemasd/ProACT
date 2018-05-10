#!/usr/bin/env python3
'''
Given a tree T, a start time S, and an end time E, remove all leaves from T that
are outside of the time window [S,E]. Trees must be in the Newick format, and
leaf labels must be delimited by the '|' character, with the last field being
the time associated with the leaf.
'''
# get the label of Dendropy node u
def label(u):
    return str(u.taxon).replace("'",'')

# run main program
if __name__ == "__main__":
    import argparse; from gzip import open as gopen; from dendropy import Tree
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
    output.write((str(tree.extract_tree_with_taxa(taxa)))); output.write('\n'); output.close()