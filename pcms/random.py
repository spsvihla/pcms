#!/usr/bin/env python3

# @file random.py
# @brief Functions related to generating random phylogenetic trees.
# @author Sean Svihla
#

import numpy as np
import ete3


## @brief Generate a random k-regular tree.
#
#  Constructs a random k-regular tree by simulating a birth and death
#  process. Note there is no guarantee that the simulation will not terminate
#  before the given size is reached.
#
#  Parameters:
#    @param int   k          The number of children descending from each node.
#    @param int   max_size   The maximum number of inernal nodes.
#    @param float lam        Birth rate.
#    @param flaot mu         Death rate.
#  
#  Outputs:
#    @return ete3.TreeNode  The phylogenetic tree.
#    @reutnr int            The number of internal nodes of the tree. 
#
def k_regular(k, max_size, lam, mu):
    rate = lam = mu
    threshold = lam / rate

    root = ete3.TreeNode()
    root.add_features(depth=0)
    node = root.add_child()
    node.add_features(depth=1)
    for __ in range(k):
        child = node.add_child(dist=np.random.exponential(rate))
        child.add_features(depth=2)
    leaves = root.get_leaves()

    n_internal_nodes = 2
    n_alive_lineages = k
    while n_internal_nodes < max_size:
        # select a random living lineage
        indx = np.random.randint(n_alive_lineages)
        leaf = leaves.pop(indx)
        n_alive_lineages += -1

        # if there is a birth, then 'leaf' becomes a parent; otherwise, 
        # the lineage is extinct
        if np.random.uniform() < threshold:
            for __ in range(k):
                child = leaf.add_child(dist=np.random.exponential(rate))
                child.add_features(depth=leaf.depth)
                leaves.append(child)
            n_alive_lineages += k
            n_internal_nodes += 1

        # terminate if all lineages go extinct
        if n_alive_lineages == 0:
            break

    return root, n_internal_nodes 
