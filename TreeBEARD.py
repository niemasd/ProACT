#!/usr/bin/env python3
# sort all (or internal) nodes by average diagnosis time and output all leaves below current node (break ties by diagnosis time)
def average(tree,diag,n,all):
    from queue import PriorityQueue,Queue; traverse = PriorityQueue()
    for u in tree.postorder_node_iter():
        u.done = False
        if u.is_leaf():
            assert u in diag, "Individual %s not in diagnostic time file" % u
            u.leaves_below = 1; u.tot_diag = diag[u]; u.avg_diag = diag[u]
        else:
            u.leaves_below = 0.; u.tot_diag = 0.
            for c in u.child_node_iter():
                u.leaves_below += c.leaves_below; u.tot_diag += c.tot_diag
            u.avg_diag = u.tot_diag/u.leaves_below
        if all or not u.is_leaf():
            traverse.put((u.avg_diag,u))
    output = []
    while not traverse.empty():
        next_avg_diag,next = traverse.get(); next.done = True
        if next.is_leaf():
            output.append(str(next.taxon).replace('"',''))
            if len(output) == n:
                return output
            continue
        to_explore = Queue(); to_explore.put((diag[next],next)); pick = PriorityQueue()
        while not to_explore.empty():
            tmp_diag,tmp = to_explore.get(); tmp.done = True
            if tmp.is_leaf():
                pick.put((tmp_diag,tmp)); continue
            for c in tmp.child_node_iter():
                if not c.done:
                    to_explore.put((diag[c],c))
        while not pick.empty():
            output.append(str(pick.get()[1].taxon).replace('"',''))
            if len(output) == n:
                return output
    assert False, "Only found %d individuals, not specified number (%d)" % (len(output),n)
def average_all(tree,diag,n):
    return average(tree,diag,n,True)
def average_internal(tree,diag,n):
    return average(tree,diag,n,False)

# run TreeBEARD
METHODS = {'average_all':average_all,'average_internal':average_internal}
if __name__ == "__main__":
    import argparse; from gzip import open as gopen; from dendropy import Tree
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=True, type=str, help="Input Tree File")
    parser.add_argument('-d', '--diagnosis', required=True, type=str, help="Input Diagnosis Time File")
    parser.add_argument('-m', '--method', required=True, type=str, help="Method")
    parser.add_argument('-n', '--number', required=True, type=int, help="Number of Individuals")
    args = parser.parse_args()
    assert args.number > 0, "Number of individuals must be a positive integer"
    assert args.method in METHODS, "Invalid method: %s. Options: %s" % (args.method, ', '.join(sorted(METHODS.keys())))
    if args.tree.lower().endswith('.gz'):
        tree = Tree.get(data=gopen(args.tree).read().decode(),schema='newick')
    else:
        tree = Tree.get(data=open(args.tree).read(),schema='newick')
    num_leaves = len(tree.leaf_nodes())
    assert args.number < num_leaves, "Number of output individuals (%d) must be less than total number of individuals in tree (%d)" % (args.number,num_leaves)
    diag = {}
    for line in {True:gopen(args.diagnosis),False:open(args.diagnosis)}[args.diagnosis.lower().endswith('.gz')]:
        if isinstance(line,bytes):
            u,t = line.decode().strip().split()
        else:
            u,t = line.strip().split()
        diag[u] = float(t)
    for u in METHODS[args.method](tree,diag,n):
        print(u)