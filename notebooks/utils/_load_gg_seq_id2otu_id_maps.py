# utils/_load_gg_seq_id2_otu_id_maps.py
import pandas as pd
from typing import Tuple


def load_gg_seq_id2otu_id_maps(filepath: str) -> Tuple[pd.Series, pd.Series]:
    """
    Load mappings from sequence ID to OTU ID and from OTU ID to list of sequence IDs.

    Parameters
    ----------
    filepath : str
        Path to the OTU map file, where each line starts with an OTU ID followed by tab-separated sequence IDs.

    Returns
    -------
    Tuple[pd.Series, pd.Series]
        - seq_id2otu_map: maps each sequence ID to its OTU ID.
        - otu2seq_id_map: maps each OTU ID to its representative sequence ID.
    """
    with open(filepath, 'r') as f:
        lines = [line.strip().split('\t') for line in f if line.strip()]

    ref_seq_ids = []
    seq_ids = []
    otu_ids = []

    for parts in lines:
        otu_ids.append(int(parts[0]))
        seq_id = list(map(int, parts[1:]))
        seq_ids.append(seq_id)
        ref_seq_ids.append(seq_id[0])

    seq_id2otu_map = pd.Series(data=otu_ids, index=ref_seq_ids)
    otu2seq_id_map = pd.Series(data=seq_ids, index=otu_ids)

    return seq_id2otu_map, otu2seq_id_map
