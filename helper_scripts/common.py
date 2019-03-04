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

def load_transmissions(filename):
    trans = []; nodes = set()
    for l in read_lines(filename):
        try:
            u,v,t = l.split(); t = float(t)
        except:
            raise RuntimeError("Invalid transmission network")
        trans.append((u,v,t))
    return trans

def load_individuals(filename):
    people = list()
    for l in read_lines(filename):
        if l.count('|') == 2: # virus|person|time identifiers
            people.append(l.split('|')[1])
        else:
            people.append(l)
    return people

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
