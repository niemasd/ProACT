#!/usr/bin/env python3
from treeswift import read_tree_newick,Node
from warnings import warn

# sort by edge length, then parent edge length, then grandparent edge length, etc.
def prioritize(tree,n):
    # compute path-length to root
    leaves = list(); root_dist = dict(); num_ancestors = dict()
    for u in tree.traverse_preorder():
        if u.is_leaf():
            leaves.append(u)
        if u.is_root():
            root_dist[u] = 0; num_ancestors[u] = 0; u.edge_length = 0
        else:
            if u.edge_length is None:
                u.edge_length = 0
            root_dist[u] = root_dist[u.parent] + u.edge_length
            num_ancestors[u] = num_ancestors[u.parent] + 1
    if n == 'All':
        n = len(leaves)
    else:
        n = int(n)
    if n < 0 or n > len(leaves):
        raise ValueError("Number of output individuals (%d) must be less than or equal to total number of individuals in tree (%d)" % (n,len(leaves)))
    # sort by edge length, then parent edge length, then grandparent edge length, etc. (use dist to root once ancestors run out)
    def compare(l1,l2):
        tie1 = float(root_dist[l1])/num_ancestors[l1]; tie2 = float(root_dist[l2])/num_ancestors[l2]
        c1 = l1; c2 = l2
        c1l = c1.edge_length; c2l = c2.edge_length
        while c1l == c2l:
            if c1 is None and c2 is None:
                return False
            if c1 is not None:
                if c1.parent is None:
                    c1l = tie1; c1 = None
                else:
                    c1 = c1.parent; c1l = c1.edge_length
            if c2 is not None:
                if c2.parent is None:
                    c2l = tie2; c2 = None
                else:
                    c2 = c2.parent; c2l = c2.edge_length
        return c1l < c2l
    Node.__lt__ = lambda self,other: compare(self,other)
    return [l.label for l in sorted(leaves)[:n]]

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
