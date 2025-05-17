# utils/_build_gg_seq_id2leaf_idx_map.py
from pcms.tree import Tree


def build_gg_seq_id2leaf_idx_map(tree: Tree) -> dict[int, int]:
    """
    Build a mapping from sequence ID (as an integer) to leaf index.

    Parameters
    ----------
    tree : Tree
        A pcms.tree.Tree instance with leaf nodes named by sequence IDs.

    Returns
    -------
    dict[int, int]
        Dictionary mapping sequence ID to its index in the leaf array.
    """
    leaves, _ = tree.find_leaves()
    return {
        int(tree.get_name(node_idx)): leaf_idx 
        for leaf_idx, node_idx in enumerate(leaves)
    }