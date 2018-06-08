#!/usr/bin/env python3
'''
Choose top n people sorted by cluster growth rate (break ties arbitrarily).
'''

# run main program
if __name__ == "__main__":
    import argparse; from gzip import open as gopen
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c', '--clustering', required=True, type=str, help="Input Clustering File (Cluster Picker format)")
    parser.add_argument('-g', '--growth', required=True, type=str, help="Input Growth Rate File")
    parser.add_argument('-n', '--number', required=True, type=int, help="Number of Individuals")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    args = parser.parse_args()
    assert args.number > 0, "Number of individuals must be a positive integer"
    if args.clustering.lower().endswith('.gz'):
        cf = gopen(args.clustering)
    else:
        cf = open(args.clustering)
    if args.growth.lower().endswith('.gz'):
        gf = gopen(args.growth)
    else:
        gf = open(args.growth)
    if args.output == 'stdout':
        from sys import stdout; output = stdout
    else:
        output = open(args.output,'w')
    cluster = {}
    for line in cf:
        if isinstance(line,bytes):
            u,c = line.decode().strip().split()
        else:
            u,c = line.strip().split()
        if u == 'SequenceName':
            continue
        cluster[u] = c
    assert args.number < len(cluster), "Number of output individuals (%d) must be less than total number of individuals" % args.number
    growth = {}
    for line in gf:
        if isinstance(line,bytes):
            c,g = line.decode().strip().split()
        else:
            c,g = line.strip().split()
        if c.startswith('Cluster') and g.startswith('GrowthRate'):
            continue
        growth[c] = float(g)
    for g,u in sorted([(growth[cluster[u]],u) for u in cluster if cluster[u] in growth], reverse=True)[:args.number]:
        output.write('%s\n'%u)
