#!/usr/bin/env python3
from common import individual_efficacy,leaf_to_name,load_diagnosis,load_transmissions
from matplotlib.cm import Reds,ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from treeswift import read_tree_newick
import argparse
import matplotlib.pyplot as plt

def compute_max_sibling_leaf_time(tree,inf):
    max_inf = dict()
    for u in tree.traverse_postorder():
        if u.is_leaf():
            max_inf[u] = inf[L2N[u]]
        else:
            max_inf[u] = max(max_inf[c] for c in u.children)
    out = dict()
    for u in tree.traverse_postorder():
        if u.is_root():
            out[u] = None
        else:
            out[u] = max(max_inf[c] for c in u.parent.children if c != u)
    return out

def edgelength_over_time(tree,inf): # each element in list for a given leaf is (time,length) tuple
    el_t = dict(); max_sib_inf = compute_max_sibling_leaf_time(tree,diag); oldest_leaf = dict()
    # compute (oldest_time,leaf) tuples for each node
    for u in tree.traverse_postorder():
        if u.is_leaf():
            oldest_leaf[u] = (inf[L2N[u]],u)
        else:
            oldest_leaf[u] = min(oldest_leaf[c] for c in u.children)
    # compute edge length over time
    for u in tree.traverse_postorder():
        if u.is_leaf():
            if inf[L2N[u]] >= max_sib_inf[u]:
                el_t[u] = [[inf[L2N[u]],u.edge_length]]
            else:
                el_t[u] = [[max_sib_inf[u],u.edge_length]]
        else:
            if u.edge_length is None:
                el = 0
            else:
                el = u.edge_length
            if u.is_root():
                prev_time = oldest_leaf[u][0]
            else:
                prev_time = min(oldest_leaf[c] for c in u.parent.children if oldest_leaf[c][1] != oldest_leaf[u][1])[0]
            if prev_time >= oldest_leaf[u][0] and prev_time < el_t[oldest_leaf[u][1]][-1][0]:
                el_t[oldest_leaf[u][1]].append([prev_time,el_t[oldest_leaf[u][1]][-1][1]])
            el_t[oldest_leaf[u][1]][-1][1] += el
    for l in el_t:
        if el_t[l][-1][0] > inf[L2N[l]]:
            el_t[l][-1][0] = inf[L2N[l]]
    return {l:el_t[l][::-1] for l in el_t}

def plot_edgelength_over_time(el_t,eff,max_num_lines):
    leaves = sorted(el_t.keys(), key=lambda a: eff[L2N[a]])
    if max_num_lines is not None:
        leaves = leaves[-max_num_lines:]
    min_eff = min(eff[L2N[l]] for l in leaves); max_eff = max(eff[L2N[l]] for l in leaves); max_time = max(el_t[l][-1][0] for l in leaves)
    norm = Normalize(vmin=min_eff, vmax=max_eff, clip=True); color_mapper = ScalarMappable(norm=norm, cmap=Reds)
    handles = [Patch(color=color_mapper.to_rgba(e), label='%d Transmission(s)'%e) for e in (min_eff,max_eff)]
    for l in leaves:
        pairs = el_t[l]
        x = [pairs[0][0]]; y = [pairs[0][1]] # start (just a point)
        for t,el in pairs[1:]:
            x.append(t); y.append(y[-1]) # bring it forward
            x.append(t); y.append(el) # bring it up
        x.append(max_time); y.append(y[-1])
        plt.plot(x, y, color=color_mapper.to_rgba(eff[L2N[l]]))
    plt.legend(handles=handles, bbox_to_anchor=(0.995, 0.995), loc=1, borderaxespad=0.)
    plt.xlabel('Time'); plt.ylabel('Edge Length'); plt.title('Edge Length vs. Time')
    plt.tight_layout(); plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tr', '--tree', required=True, type=str, help="Phylogenetic Tree (Newick format)")
    parser.add_argument('-d', '--diagnosis', required=True, type=str, help="Diagnosis Times (TSV)")
    parser.add_argument('-tn', '--transmissions', required=True, type=str, help="Transmission Network (FAVITES format)")
    parser.add_argument('-t', '--from_time', required=True, type=float, help="From Time (for # transmissions)")
    parser.add_argument('-tt', '--to_time', required=False, type=float, default=float('inf'), help="To Time (for # transmissions)")
    parser.add_argument('-n', '--max_num_lines', required=False, type=int, default=None, help="Maximum Number of Lines to Draw (reduce this for speed)")
    args = parser.parse_args()
    tree = read_tree_newick(args.tree)
    diag = load_diagnosis(args.diagnosis)
    global L2N; L2N = leaf_to_name(tree)
    eff = individual_efficacy([L2N[l] for l in tree.traverse_leaves()],load_transmissions(args.transmissions),args.from_time,args.to_time)
    el_t = edgelength_over_time(tree,diag)
    plot_edgelength_over_time(el_t,eff,args.max_num_lines)
