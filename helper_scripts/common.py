#!/usr/bin/env python3
'''
Common functions
'''
from gzip import open as gopen
from sys import stdin

def read_lines(filename):
    if filename == 'stdin':
        return [l.strip() for l in stdin.read().strip().splitlines()]
    elif filename.lower().endswith('.gz'):
        return [l.strip() for l in gopen(filename).read().decode().strip().splitlines()]
    else:
        return [l.strip() for l in open(filename).read().strip().splitlines()]

def load_diag_times(filename):
    diag = dict()
    for l in read_lines(filename):
        u,t = l.strip().split('\t')
        if u.count('|') == 2: # virus|person|time identifiers
            u = u.split('|')[1]
        diag[u] = float(t)
    return diag

def load_transmissions(filename):
    trans = []; nodes = set()
    for l in read_lines(filename):
        try:
            u,v,t = l.split(); t = float(t)
        except:
            raise RuntimeError("Invalid transmission network")
        trans.append((u,v,t))
    return trans

def individuals_from_lines(lines):
    people = list()
    for l in lines:
        if l.count('|') == 2: # virus|person|time identifiers
            people.append(l.split('|')[1])
        else:
            people.append(l)
    return people

def leaf_to_name(tree):
    tr = dict()
    for l in tree.traverse_leaves():
        if l.label.count('|') == 2: # virus|person|time identifiers
            tr[l] = l.label.split('|')[1]
        else:
            tr[l] = l.label
    return tr

def load_individuals(filename):
    return individuals_from_lines(read_lines(filename))

def individuals_from_tree(tree):
    return individuals_from_lines([l.label for l in tree.traverse_leaves()])

def individual_efficacy(user_individuals,transmissions,from_time,to_time):
    assert to_time > from_time, "To Time must be larger than From Time"
    eff = {u:0 for u in user_individuals}; trans_nodes = set()
    for u,v,t in transmissions:
        trans_nodes.add(u); trans_nodes.add(v)
        if t >= from_time and t <= to_time and u in eff:
            eff[u] += 1
    for u in eff:
        assert u in trans_nodes, "Individual not in transmission network: %s"%u
    return eff

def optimal_order(individuals,eff):
    return sorted(individuals, key=lambda x: eff[x], reverse=True)
