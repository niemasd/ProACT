#!/usr/bin/env python3
from dendropy import Node,Tree
from queue import PriorityQueue,Queue
from warnings import warn
UNUSED_DIAG_WARNING = "Specified diagnostic time file, but it will be ignored in this method"

# get the label of Dendropy node u
def label(u):
    return str(u.taxon).replace("'",'')

# create < operator for Dendropy Node
def node_lt(self, other):
    return label(self) < label(other)
Node.__lt__ = node_lt

# randomly pick n individuals
def random_select(tree,diag,n):
    if diag is not None:
        warn(UNUSED_DIAG_WARNING)
    from random import shuffle
    leaves = [label(leaf) for leaf in tree.leaf_node_iter()]
    shuffle(leaves)
    return leaves[:n]

# sort all (or internal) nodes by average diagnosis time and output all leaves below current node (break ties by diagnosis time)
def average(tree,diag,n,all):
    traverse = PriorityQueue()
    for u in tree.postorder_node_iter():
        u.done = False
        if u.is_leaf():
            u_label = label(u)
            assert u_label in diag, "Individual %s not in diagnostic time file" % u_label
            u.leaves_below = 1; u.tot_diag = diag[u_label]; u.avg_diag = diag[u_label]
        else:
            u.leaves_below = 0.; u.tot_diag = 0.
            for c in u.child_node_iter():
                u.leaves_below += c.leaves_below; u.tot_diag += c.tot_diag
            u.avg_diag = u.tot_diag/u.leaves_below
        if all or not u.is_leaf():
            traverse.put((-u.avg_diag,u))
    output = []
    while not traverse.empty():
        dummy,next = traverse.get(); next.done = True
        if next.is_leaf():
            output.append(label(next))
            if len(output) == n:
                return output
            continue
        to_explore = Queue(); to_explore.put(next); pick = PriorityQueue()
        while not to_explore.empty():
            tmp = to_explore.get(); tmp.done = True
            if tmp.is_leaf():
                pick.put((-diag[label(tmp)],tmp)); continue
            for c in tmp.child_node_iter():
                if not c.done:
                    to_explore.put(c)
        while not pick.empty():
            output.append(label(pick.get()[1]))
            if len(output) == n:
                return output
    assert False, "Only found %d individuals, not specified number (%d)" % (len(output),n)
def average_all(tree,diag,n):
    return average(tree,diag,n,True)
def average_internal(tree,diag,n):
    return average(tree,diag,n,False)

# run TreeBEARD
METHODS = {'average_all':average_all,'average_internal':average_internal,'random':random_select}
NEED_DIAGNOSIS = {'average_all','average_internal'}
if __name__ == "__main__":
    import argparse; from gzip import open as gopen
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, help="Input Tree File")
    parser.add_argument('-d', '--diagnosis', required=False, type=str, default=None, help="Input Diagnosis Time File")
    parser.add_argument('-m', '--method', required=True, type=str, help="Method")
    parser.add_argument('-n', '--number', required=True, type=int, help="Number of Individuals")
    args = parser.parse_args()
    assert args.number > 0, "Number of individuals must be a positive integer"
    args.method = args.method.lower()
    assert args.method in METHODS, "Invalid method: %s. Options: %s" % (args.method, ', '.join(sorted(METHODS.keys())))
    if args.diagnosis is None:
        assert args.method not in NEED_DIAGNOSIS, "Method %s requires diagnosis file" % args.method
        diag = None
    else:
        diag = {}
        for line in {True:gopen(args.diagnosis),False:open(args.diagnosis)}[args.diagnosis.lower().endswith('.gz')]:
            if isinstance(line,bytes):
                u,t = line.decode().strip().split()
            else:
                u,t = line.strip().split()
            diag[u] = float(t)
    if args.tree.lower().endswith('.gz'):
        tree = Tree.get(data=gopen(args.tree).read().decode(),schema='newick')
    else:
        tree = Tree.get(data=open(args.tree).read(),schema='newick')
    num_leaves = len(tree.leaf_nodes())
    assert args.number < num_leaves, "Number of output individuals (%d) must be less than total number of individuals in tree (%d)" % (args.number,num_leaves)
    for u in METHODS[args.method](tree,diag,args.number):
        print(u)