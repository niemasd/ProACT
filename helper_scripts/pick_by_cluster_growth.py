#!/usr/bin/env python3
'''
Choose top n people sorted by cluster growth rate (break ties arbitrarily).
'''

# run main program
if __name__ == "__main__":
    import argparse; from gzip import open as gopen; from warnings import warn
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c', '--clustering', required=True, type=str, help="Input Clustering File (Cluster Picker format)")
    parser.add_argument('-g', '--growth', required=True, type=str, help="Input Growth Rate File")
    parser.add_argument('-d', '--diagnosis', required=False, type=str, default=None, help="Diagnosis File")
    parser.add_argument('-n', '--number', required=False, type=str, default='All', help="Number of Individuals")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    args = parser.parse_args()
    if args.clustering.lower().endswith('.gz'):
        cf = gopen(args.clustering)
    else:
        cf = open(args.clustering)
    if args.growth.lower().endswith('.gz'):
        gf = gopen(args.growth)
    else:
        gf = open(args.growth)
    if args.diagnosis is not None:
        if args.diagnosis.lower().endswith('.gz'):
            df = gopen(args.diagnosis).read().decode().strip().splitlines()
        else:
            df = open(args.diagnosis).read().strip().splitlines()
        everybody = {l.split('\t')[0].strip() for l in df}
    else:
        everybody = set()
    if args.output == 'stdout':
        from sys import stdout; output = stdout
    else:
        output = open(args.output,'w')
    cluster = {}; num_people = 0
    for line in cf:
        if isinstance(line,bytes):
            u,c = line.decode().strip().split()
        else:
            u,c = line.strip().split()
        if u == 'SequenceName':
            continue
        cluster[u] = c; num_people += 1; everybody.add(u)
    if args.number == 'All':
        args.number = len(everybody)
    else:
        args.number = int(args.number)
    assert args.number > 0, "Number of individuals must be a positive integer"
    if args.number > len(cluster):
        assert args.diagnosis is not None, "Number of output individuals (%d) is greater than the total number of individuals (%d), so must specify diagnosis file" % (args.number, len(cluster))
        warn("Number of output individuals (%d) is greater than the total number of individuals (%d), so the remaining %d individuals will be randomly selected from the diagnosis file (%s)." % (args.number, len(cluster), args.number-len(cluster), args.diagnosis))
    growth = {}
    for line in gf:
        if isinstance(line,bytes):
            c,g = line.decode().strip().split()
        else:
            c,g = line.strip().split()
        if c.startswith('Cluster') and g.startswith('GrowthRate'):
            continue
        growth[c] = float(g)
    num_output = 0
    for g,u in sorted([(growth[cluster[u]],u) for u in cluster if cluster[u] in growth], reverse=True)[:args.number]:
        output.write('%s\n'%u); everybody.remove(u); num_output += 1
    while num_output < args.number:
        output.write('%s\n'%everybody.pop()); num_output += 1
