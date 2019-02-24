#!/usr/bin/env python3
from treeswift import read_tree_newick,Node
from random import shuffle
from queue import PriorityQueue,Queue
from warnings import warn
import heapq
USE_TREE_WARNING = "Only tree file will be used"
USE_INF_WARNING = "Only diagnosis time file will be used"
NEED_TREE_INF_ERROR = "User did not specify tree nor diagnosis time file, but one of the two is needed in this method"
UNUSED_TREE_WARNING = "User specified tree file, but it will be ignored in this method"
NEED_TREE_ERROR = "No tree file specified, but it is needed for this method"
UNUSED_INF_WARNING = "User specified diagnosis time file, but it will be ignored in this method"
NEED_INF_ERROR = "No diagnosis time file specified, but it is needed for this method"
INVALID_DATE = "Invalid date. Dates must be floats/integers"

# randomly pick n individuals
def random_select(tree,inf,n):
    assert tree is not None or inf is not None, NEED_TREE_INF_ERROR
    if tree is not None and inf is not None:
        warn(USE_TREE_WARNING)
    if tree is not None:
        individuals = [str(leaf) for leaf in tree.traverse_leaves()]
    elif inf is not None:
        individuals = list(inf.keys())
    shuffle(individuals)
    return individuals[:n]

# feedback mode
def feedback(tree,inf,n):
    assert tree is not None, NEED_TREE_ERROR
    assert inf is not None, NEED_INF_ERROR
    pq = list(); output = list(); outset = set()
    for u in tree.traverse_postorder():
        u.done = False
        if u.is_leaf():
            if str(u) not in inf:
                warn("Individual %s not in diagnosis time file, so assuming diagnosis time = infinity" % str(u))
                inf[str(u)] = float('inf')
            u.num = inf[str(u)]; u.den = 1.
        else:
            u.num = 0.; u.den = 0.
            for c in u.children:
                u.num += c.num; u.den += c.den
            u.tup = [-1*u.num/u.den, u]; pq.append(u.tup)
    heapq.heapify(pq)
    while len(pq) != 0 and len(output) < n:
        u = heapq.heappop(pq)[1]
        if u.done:
            continue
        leaves_below = list()
        for d in u.traverse_preorder():
            d.done = True
            if d.is_leaf():
                leaves_below.append((-1*inf[str(d)],d))
        leaves_below.sort()
        to_fix = set()
        for DUMMY,u in leaves_below:
            ustr = str(u)
            if ustr not in outset:
                output.append(ustr); outset.add(ustr)
            for a in u.traverse_ancestors(include_self=False):
                to_fix.add(u); a.num -= u.num; a.den -= 1.
        for u in to_fix:
            if not u.is_leaf():
                u.tup[0] = -1*u.num/u.den
        heapq.heapify(pq)
    return output[:n]

# sort all (or internal) nodes by average diagnosis time and output all leaves below current node (break ties by diagnosis time)
def average(tree,inf,n,sort_max,all):
    assert tree is not None, NEED_TREE_ERROR
    assert inf is not None, NEED_INF_ERROR
    traverse = PriorityQueue()
    for u in tree.traverse_postorder():
        u.done = False
        if u.is_leaf():
            if str(u) not in inf:
                warn("Individual %s not in diagnosis time file, so assuming diagnosis time = infinity" % str(u))
                inf[str(u)] = float('inf')
            u.leaves_below = 1; u.tot_inf = inf[str(u)]; u.avg_inf = inf[str(u)]
        else:
            u.leaves_below = 0.; u.tot_inf = 0.
            for c in u.children:
                u.leaves_below += c.leaves_below; u.tot_inf += c.tot_inf
            u.avg_inf = u.tot_inf/u.leaves_below
        if all or not u.is_leaf():
            traverse.put(({True:-u.avg_inf,False:u.avg_inf}[sort_max],u))
    output = list(); outset = set()
    while not traverse.empty():
        dummy,next = traverse.get(); next.done = True
        if next.is_leaf():
            nextstr = str(next)
            if nextstr not in outset:
                output.append(nextstr); outset.add(nextstr)
            if len(output) == n:
                return output
            continue
        to_explore = Queue(); to_explore.put(next); pick = PriorityQueue()
        while not to_explore.empty():
            tmp = to_explore.get(); tmp.done = True
            if tmp.is_leaf():
                pick.put(({True:-inf[str(tmp)],False:inf[str(tmp)]}[sort_max],tmp)); continue
            for c in tmp.children:
                if not c.done:
                    to_explore.put(c)
        while not pick.empty():
            ustr = str(pick.get()[1])
            if ustr not in outset:
                output.append(ustr); outset.add(ustr)
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

# sort all individuals by diagnosis time
def sort_by_inf(tree,inf,n,sort_max):
    if tree is not None:
        warn(UNUSED_TREE_WARNING)
    return [e[0] for e in sorted(inf.items(), key=lambda k: (k[1],k[0]), reverse=sort_max)[:n]]
def max_inf(tree,inf,n):
    return sort_by_inf(tree,inf,n,True)
def min_inf(tree,inf,n):
    return sort_by_inf(tree,inf,n,False)

# sort all individuals by distance from root
def sort_by_root_dist(tree,inf,n,sort_max):
    if inf is not None:
        warn(UNUSED_INF_WARNING)
    return [str(l[1]) for l in tree.traverse_rootdistorder(ascending=sort_max) if l[1].is_leaf()][:n]
def max_root_dist(tree,inf,n):
    return sort_by_root_dist(tree,inf,n,False)
def min_root_dist(tree,inf,n):
    return sort_by_root_dist(tree,inf,n,True)

# run ProACT
METHODS = {
    'average_max_inf_all':average_max_inf_all,
    'average_max_inf_internal':average_max_inf_internal,
    'average_min_inf_all':average_min_inf_all,
    'average_min_inf_internal':average_min_inf_internal,
    'feedback':feedback,
    'max_inf':max_inf,
    'min_inf':min_inf,
    'max_root_dist':max_root_dist,
    'min_root_dist':min_root_dist,
    'random':random_select}
if __name__ == "__main__":
    import argparse; from gzip import open as gopen
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str, default=None, help="Input Tree File")
    parser.add_argument('-d', '--diagnosis', required=False, type=str, default=None, help="Input Diagnosis Time File")
    parser.add_argument('-m', '--method', required=False, type=str, default='average_max_inf_internal', help="Method (%s)" % ', '.join(sorted(METHODS.keys())))
    parser.add_argument('-n', '--number', required=False, type=str, default='All', help="Number of Individuals")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    args = parser.parse_args()
    args.method = args.method.lower()
    assert args.method in METHODS, "Invalid method: %s. Options: %s" % (args.method, ', '.join(sorted(METHODS.keys())))
    assert args.diagnosis is not None or args.tree is not None, "Must specify either tree file or diagnosis time file (or both)"
    if args.output == 'stdout':
        from sys import stdout; output = stdout
    else:
        output = open(args.output,'w')
    if args.diagnosis is None:
        inf = None
    else:
        inf = {}
        for line in {True:gopen(args.diagnosis),False:open(args.diagnosis)}[args.diagnosis.lower().endswith('.gz')]:
            if isinstance(line,bytes):
                u,t = line.decode().strip().split()
            else:
                u,t = line.strip().split()
            try:
                inf[u] = float(t)
            except ValueError:
                raise ValueError(INVALID_DATE)
    if args.tree is None:
        tree = None
        if args.number == 'All':
            args.number = len(inf)
        else:
            args.number = int(args.number)
    else:
        if args.tree.lower().endswith('.gz'):
            tree = read_tree_newick(gopen(args.tree).read().decode())
        else:
            tree = read_tree_newick(open(args.tree).read())
        num_leaves = len([l for l in tree.traverse_leaves()])
        if args.number == 'All':
            args.number = num_leaves
        else:
            args.number = int(args.number)
        assert 0 < args.number <= num_leaves, "Number of output individuals (%d) must be less than or equal to total number of individuals in tree (%d)" % (args.number,num_leaves)
    for u in METHODS[args.method](tree,inf,args.number):
        output.write(u); output.write('\n')
