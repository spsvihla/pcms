# utils/_load_gg_seq_id2_otu_id_maps.py
from typing import Tuple


def load_gg_seq_id2otu_id_maps(filepath: str) -> Tuple[dict[int, int], dict[int, list[int]]]:
    """
    Load mappings from sequence ID to OTU ID and from OTU ID to list of sequence IDs.

    Parameters
    ----------
    filepath : str
        Path to the OTU map file, where each line starts with an OTU ID followed by tab-separated sequence IDs.

    Returns
    -------
    Tuple[dict[int, int], dict[int, list[int]]]
        - seq_id2otu_map: maps each sequence ID to its OTU ID.
        - otu2seq_id_map: maps each OTU ID to its representative sequence ID.
    """
    with open(filepath, 'r') as f:
        lines = [line.strip().split('\t') for line in f if line.strip()]

    seq_id2otu_map = {}
    otu2seq_id_map = {}

    for parts in lines:
        otu_id = int(parts[0])
        seq_ids = list(map(int, parts[1:]))
        otu2seq_id_map[otu_id] = seq_ids[0]
        seq_id2otu_map.update({seq_id: otu_id for seq_id in seq_ids})

    return seq_id2otu_map, otu2seq_id_map
