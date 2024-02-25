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
#    @param int   max_size   The maximum number of internal nodes.
#    @param float lam        Birth rate.
#    @param flaot mu         Death rate.
#  
#  Outputs:
#    @return ete3.TreeNode  The phylogenetic tree.
#    @return int            The number of leaves in the tree.
#
def k_regular(k, max_size, lam, mu):
    rate = lam = mu
    threshold = lam / rate

    root = ete3.TreeNode()
    node = root.add_child()
    for __ in range(k):
        node.add_child(dist=np.random.exponential(rate))
    leaves = root.get_leaves()

    size = 1
    n_leaves = k
    while size < max_size:
        # select a random living lineage
        indx = np.random.randint(n_leaves)
        leaf = leaves.pop(indx)
        n_leaves += -1

        # if there is a birth, then 'leaf' becomes a parent; otherwise, 
        # the lineage is extinct
        if np.random.uniform() < threshold:
            for __ in range(k):
                leaves.append(leaf.add_child(dist=np.random.exponential(rate)))
            n_leaves += k
            size += 1

        # terminate if all lineages go extinct
        if n_leaves == 0:
            break

    return root
