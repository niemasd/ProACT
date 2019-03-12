#!/usr/bin/env python3
from treeswift import read_tree_newick,Node
from warnings import warn

# sort by edge length, then parent edge length, then grandparent edge length, etc.
def prioritize(tree,n):
    # compute path-length to root
    leaves = list(); root_dist = dict()#; root_dist_u = dict()
    for u in tree.traverse_preorder():
        if u.is_root():
            root_dist[u] = 0#; root_dist_u[u] = 0
        else:
            if u.edge_length is None:
                u.edge_length = 0
            root_dist[u] = root_dist[u.parent] + u.edge_length#; root_dist_u[u] = root_dist_u[u.parent] + 1
        if u.is_leaf():
            leaves.append(u)
    avg_root_dist = sum(root_dist[u] for u in leaves)/len(leaves)
    if n == 'All':
        n = len(leaves)
    else:
        n = int(n)
    if n < 0 or n > len(leaves):
        raise ValueError("Number of output individuals (%d) must be less than or equal to total number of individuals in tree (%d)" % (n,len(leaves)))
    return [l.label for l in sorted(leaves, key=lambda x: (x.edge_length, abs(root_dist[x]-avg_root_dist)))[:n]]

# run ProACT
if __name__ == "__main__":
    import argparse; from gzip import open as gopen
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, help="Input Tree File")
    parser.add_argument('-n', '--number', required=False, type=str, default='All', help="Number of Individuals")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    args = parser.parse_args()
    if args.output == 'stdout':
        from sys import stdout; output = stdout
    else:
        output = open(args.output,'w')
    tree = read_tree_newick(args.tree)
    for u in prioritize(tree,args.number):
        output.write(u); output.write('\n')
