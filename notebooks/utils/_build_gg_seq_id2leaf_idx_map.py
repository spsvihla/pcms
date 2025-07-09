# utils/_build_gg_seq_id2leaf_idx_map.py
import pandas as pd
from pcms.tree import Tree


def build_gg_seq_id2leaf_idx_map(tree: Tree) -> pd.Series:
    """
    Build a mapping from sequence ID (as an integer) to leaf index.

    Parameters
    ----------
    tree : Tree
        A pcms.tree.Tree instance with leaf nodes named by sequence IDs.

    Returns
    -------
    pd.Series
        Leaf indexes indexed by sequence id.
    """
    leaves, _ = tree.find_leaves()
    seq_ids = []
    leaf_idxs = []
    for leaf_idx, node_idx in enumerate(leaves):
        seq_ids.append(tree.get_name(node_idx))
        leaf_idxs.append(leaf_idx)
    seq_id2leaf_idx_map = pd.Series(data=leaf_idxs, index=seq_ids)
    return seq_id2leaf_idx_map