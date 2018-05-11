#!/usr/bin/env python3
from dendropy import Node,Tree
from random import shuffle
from queue import PriorityQueue,Queue
from warnings import warn
TREE_OVER_INF_WARNING = "User specified both tree and infection time files, but only tree file will be used"
NEED_TREE_OR_INF_ERROR = "User did not specify tree nor infection time file, but one of the two is needed in this method"
UNUSED_TREE_WARNING = "User specified tree file, but it will be ignored in this method"
NEED_TREE_ERROR = "No tree file specified, but it is needed for this method"
UNUSED_INF_WARNING = "User specified infection time file, but it will be ignored in this method"
NEED_INF_ERROR = "No infection time file specified, but it is needed for this method"

# get the label of Dendropy node u
def label(u):
    return str(u.taxon).replace("'",'')

# create < operator for Dendropy Node
def node_lt(self, other):
    return label(self) < label(other)
Node.__lt__ = node_lt

# randomly pick n individuals
def random_select(tree,inf,n):
    assert tree is not None or inf is not None, NEED_TREE_OR_INF_ERROR
    if tree is not None and inf is not None:
        warn(TREE_OVER_INF_WARNING)
    if tree is not None:
        individuals = [label(leaf) for leaf in tree.leaf_node_iter()]
    elif inf is not None:
        individuals = list(inf.keys())
    shuffle(individuals)
    return individuals[:n]


# sort all (or internal) nodes by average infection time and output all leaves below current node (break ties by infection time)
def average(tree,inf,n,sort_max,all):
    assert tree is not None, NEED_TREE_ERROR
    assert inf is not None, NEED_INF_ERROR
    traverse = PriorityQueue()
    for u in tree.postorder_node_iter():
        u.done = False
        if u.is_leaf():
            u_label = label(u)
            assert u_label in inf, "Individual %s not in infection time file" % u_label
            u.leaves_below = 1; u.tot_inf = inf[u_label]; u.avg_inf = inf[u_label]
        else:
            u.leaves_below = 0.; u.tot_inf = 0.
            for c in u.child_node_iter():
                u.leaves_below += c.leaves_below; u.tot_inf += c.tot_inf
            u.avg_inf = u.tot_inf/u.leaves_below
        if all or not u.is_leaf():
            traverse.put(({True:-u.avg_inf,False:u.avg_inf}[sort_max],u))
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
                pick.put(({True:-inf[label(tmp)],False:inf[label(tmp)]}[sort_max],tmp)); continue
            for c in tmp.child_node_iter():
                if not c.done:
                    to_explore.put(c)
        while not pick.empty():
            output.append(label(pick.get()[1]))
            if len(output) == n:
                return output
    assert False, "Only found %d individuals, not specified number (%d)" % (len(output),n)
def average_max_inf_all(tree,inf,n):
    return average(tree,inf,n,True,True)
def average_max_inf_internal(tree,inf,n):
    return average(tree,inf,n,True,False)
def average_min_inf_all(tree,inf,n):
    return average(tree,inf,n,False,True)
def average_min_inf_internal(tree,inf,n):
    return average(tree,inf,n,False,False)

# sort all individuals by infection time
def sort_by_inf(tree,inf,n,sort_max):
    if tree is not None:
        warn(UNUSED_TREE_WARNING)
    return [e[0] for e in sorted(inf.items(), key=lambda k: (k[1],k[0]), reverse=sort_max)[:n]]
def max_inf(tree,inf,n):
    return sort_by_inf(tree,inf,n,True)
def min_inf(tree,inf,n):
    return sort_by_inf(tree,inf,n,False)

# run TreeBEARD
METHODS = {
    'average_max_inf_all':average_max_inf_all,
    'average_max_inf_internal':average_max_inf_internal,
    'average_min_inf_all':average_min_inf_all,
    'average_min_inf_internal':average_min_inf_internal,
    'max_inf':max_inf,
    'min_inf':min_inf,
    'random':random_select}
if __name__ == "__main__":
    import argparse; from gzip import open as gopen
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str, default=None, help="Input Tree File")
    parser.add_argument('-i', '--infection', required=False, type=str, default=None, help="Input Infection Time File")
    parser.add_argument('-m', '--method', required=True, type=str, help="Method (%s)" % ', '.join(sorted(METHODS.keys())))
    parser.add_argument('-n', '--number', required=True, type=int, help="Number of Individuals")
    args = parser.parse_args()
    assert args.number > 0, "Number of individuals must be a positive integer"
    args.method = args.method.lower()
    assert args.method in METHODS, "Invalid method: %s. Options: %s" % (args.method, ', '.join(sorted(METHODS.keys())))
    if args.infection is None:
        inf = None
    else:
        inf = {}
        for line in {True:gopen(args.infection),False:open(args.infection)}[args.infection.lower().endswith('.gz')]:
            if isinstance(line,bytes):
                u,t = line.decode().strip().split()
            else:
                u,t = line.strip().split()
            inf[u] = float(t)
    if args.tree is None:
        tree = None
    else:
        if args.tree.lower().endswith('.gz'):
            tree = Tree.get(data=gopen(args.tree).read().decode(),schema='newick')
        else:
            tree = Tree.get(data=open(args.tree).read(),schema='newick')
        num_leaves = len(tree.leaf_nodes())
        assert args.number < num_leaves, "Number of output individuals (%d) must be less than total number of individuals in tree (%d)" % (args.number,num_leaves)
    for u in METHODS[args.method](tree,inf,args.number):
        print(u)