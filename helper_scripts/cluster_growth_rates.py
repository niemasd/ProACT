#!/usr/bin/env python3
'''
Given two clusterings in the Cluster Picker format, output cluster growth rates.
'''

# growth rate function
def growth(n1,n2):
    return (n2-n1)/(n2**0.5)

# run main program
if __name__ == "__main__":
    import argparse; from gzip import open as gopen
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c1', '--clustering1', required=True, type=str, help="Input Clustering 1 File (Cluster Picker format)")
    parser.add_argument('-c2', '--clustering2', required=True, type=str, help="Input Clustering 2 File (Cluster Picker format)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File")
    args = parser.parse_args()
    if args.clustering1.lower().endswith('.gz'):
        c1 = gopen(args.clustering1)
    else:
        c1 = open(args.clustering1)
    if args.clustering2.lower().endswith('.gz'):
        c2 = gopen(args.clustering2)
    else:
        c2 = open(args.clustering2)
    if args.output == 'stdout':
        from sys import stdout; output = stdout
    else:
        output = open(args.output,'w')
    n = {'c1':{},'c2':{}}
    for f,k in [(c1,'c1'),(c2,'c2')]:
        for line in f:
            if isinstance(line,bytes):
                u,c = line.decode().strip().split()
            else:
                u,c = line.strip().split()
            if u == 'SequenceName':
                continue
            if c not in n[k]:
                n[k][c] = 0
            n[k][c] += 1
    output.write('Cluster\tGrowthRate\n')
    for c in n['c2']:
        if c not in n['c1']:
            n['c1'][c] = 0
        output.write('%s\t%f\n'%(c,growth(n['c1'][c],n['c2'][c])))
