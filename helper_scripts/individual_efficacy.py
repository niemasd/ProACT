#!/usr/bin/env python3
'''
Compute clustering efficacy (average number of individuals infected by
user-selected individuals between from_time and to_time).
'''
from common import individual_efficacy,load_individuals,load_transmissions
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--individuals', required=False, type=str, default='stdin', help="Individuals (one per line)")
    parser.add_argument('-tn', '--transmissions', required=True, type=str, help="Transmission Network (FAVITES format)")
    parser.add_argument('-t', '--from_time', required=True, type=float, help="From Time")
    parser.add_argument('-tt', '--to_time', required=False, type=float, default=float('inf'), help="To Time")
    args = parser.parse_args()
    trans = load_transmissions(args.transmissions)
    user_individuals = load_individuals(args.individuals)
    eff = individual_efficacy(user_individuals,trans,args.from_time,args.to_time)
    for u in user_individuals:
        print('%s\t%d'%(u,eff[u]))
