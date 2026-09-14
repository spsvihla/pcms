from typing import Optional, Tuple
from numpy.typing import NDArray
import numpy as np
import pandas as pd
from pcms.tree import Tree


##
# build_gg_otu_id2tax_map
##
def _load_tax_filepath(tax_filepath: str) -> pd.DataFrame:
    try:
        table = pd.read_table(tax_filepath, sep='\t', names=['OTU ID', 'Taxonomy'], index_col=0)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Missing file: {e.filename}") from e
    except pd.errors.ParserError as e:
        raise RuntimeError(f"Parsing error: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Unexpected error while loading files: {e}") from e
    return table


def _get_lowest_common_rank(taxa: pd.DataFrame) -> str:
    tax_split = [row['Taxonomy'].split(';') for _, row in taxa.iterrows()]
    ranks_by_level = list(zip(*tax_split))
    lca = []
    for level in ranks_by_level:
        if len(set(level)) == 1 and level[0]:
            lca.append(level[0])
        else:
            break
    return ';'.join(lca)


def _get_relative_abundances_below_lcr(
    taxa: pd.DataFrame, lcr: str, a_abunds: NDArray, b_abunds: NDArray = None
) -> pd.DataFrame:
    if not lcr:
        return pd.DataFrame(columns=['otu_fracs', 'abund_fracs'], dtype=float)

    lcr_levels = lcr.split(';')
    lcr_depth = len(lcr_levels)

    # determine taxa below lowest common rank
    next_level_taxa = []
    for tax in taxa['Taxonomy']:
        split_tax = tax.split(';')
        if len(split_tax) > lcr_depth:
            next_level_taxa.append(split_tax[lcr_depth])
        else:
            next_level_taxa.append('unclassified')

    # OTU-count fractions
    otu_fracs = pd.Series(next_level_taxa).value_counts(normalize=True)

    # abundance-weighted fractions
    df = pd.DataFrame({'taxon': next_level_taxa, 'a_abund': a_abunds, 'b_abund': b_abunds})
    a_abund_fracs = df.groupby('taxon')['a_abund'].sum()

    if b_abunds is None:
        result = pd.concat([otu_fracs, a_abund_fracs], axis=1)
        result.columns = ['otu_fracs', 'abund_fracs']
    else:
        b_abund_fracs = df.groupby('taxon')['b_abund'].sum()
        result = pd.concat([otu_fracs, a_abund_fracs, b_abund_fracs], axis=1)
        result.columns = ['otu_fracs', 'a_abund_fracs', 'b_abund_fracs']

    return result


def build_gg_otu_id2tax_map(
    tree: Tree, tax_filepath: str, node: int, a_abunds: NDArray, b_abunds: Optional[NDArray] = None
) -> dict:
    table = _load_tax_filepath(tax_filepath=tax_filepath)
    leaves = tree.find_leaves()
    otus = np.array([tree.get_name(i) for i in leaves]).astype(int)
    start = tree.find_subtree_start_indices()[node]
    size = tree.get_subtree_size()[node]

    taxa = table.loc[otus[start:start+size]]
    lcr = _get_lowest_common_rank(taxa=taxa)
    if b_abunds is None:
        fracs = _get_relative_abundances_below_lcr(taxa, lcr, a_abunds[start:start+size])
    else:
        fracs = _get_relative_abundances_below_lcr(taxa, lcr, a_abunds[start:start+size], b_abunds[start:start+size])
    
    return {'lcr': lcr, 'fracs': fracs}


##
# build_gg_seq_id2leaf_idx_map
##
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
    leaves = tree.find_leaves()
    seq_ids = []
    leaf_idxs = []
    for leaf_idx, node_idx in enumerate(leaves):
        seq_ids.append(int(tree.get_name(node_idx)))
        leaf_idxs.append(leaf_idx)
    seq_id2leaf_idx_map = pd.Series(data=leaf_idxs, index=seq_ids)
    return seq_id2leaf_idx_map


##
# load_bs_samp2site_map
##
def load_bs_samp2site_map(
    table_filepath: str
) -> pd.DataFrame:
    """
    Load sample-to-body-site map.

    Parameters
    ----------
    table_filepath : str
        Path to the TSV file containing mapping.

    Returns
    -------
    pd.DataFrame
        DataFrame of mapping.
    """
    try:
        data = pd.read_table(table_filepath, sep='\t', engine='c')
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Missing file: {e.filename}") from e
    except pd.errors.ParserError as e:
        raise RuntimeError(f"Parsing error: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Unexpected error while loading files: {e}") from e

    data.drop('No.', axis=1)

    grouped = data.groupby('Body habitat')['SampleID'].apply(lambda x: [s.upper() for s in x])
    return grouped


##
# load_gg_seq_id2otu_id_maps
##
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


##
# load_seq_id_abundances_from_biom_tsv
##
def load_seq_id_abundances_from_biom_tsv(
    table_filepath: str,
) -> pd.DataFrame:
    """
    Load sequence abundance data and associated metadata from TSV files.

    Parameters
    ----------
    table_filepath : str
        Path to the TSV file containing sequence abundance data.

    Returns
    -------
    pd.DataFrame
        DataFrame of abundance values (rows are sequence IDs, columns are samples)
    """
    try:
        data = pd.read_table(table_filepath, skiprows=1, header=0, sep='\t', engine='c')
        data.rename(columns={'#OTU ID': 'OTU ID'}, inplace=True)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Missing file: {e.filename}") from e
    except pd.errors.ParserError as e:
        raise RuntimeError(f"Parsing error: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Unexpected error while loading files: {e}") from e

    return data