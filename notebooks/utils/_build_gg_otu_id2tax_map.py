import numpy as np
import pandas as pd
from pcms.tree import Tree


def load_tax_filepath(tax_filepath: str) -> pd.DataFrame:
    try:
        table = pd.read_table(tax_filepath, sep='\t', names=['OTU ID', 'Taxonomy'], index_col=0)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Missing file: {e.filename}") from e
    except pd.errors.ParserError as e:
        raise RuntimeError(f"Parsing error: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Unexpected error while loading files: {e}") from e
    return table


def get_lowest_common_rank(taxonomies: pd.DataFrame) -> str:
    tax_split = [row['Taxonomy'].split(';') for _, row in taxonomies.iterrows()]
    ranks_by_level = list(zip(*tax_split))
    lca = []
    for level in ranks_by_level:
        if len(set(level)) == 1 and level[0]:
            lca.append(level[0])
        else:
            break
    return ';'.join(lca)


def get_relative_abundances_below_lcr(
    taxonomies: pd.DataFrame, lcr: str, abunds: np.ndarray
) -> pd.DataFrame:
    """
    Given a set of taxonomy strings and the LCR, compute both the reference
    (count-based) fractions and the abundance-weighted fractions of the taxa
    at the rank just below the LCR.

    Parameters
    ----------
    taxonomies : pd.DataFrame
        DataFrame with a 'Taxonomy' column containing semicolon-separated taxonomy strings.
    lcr : str
        Lowest common rank string.
    abunds : np.ndarray
        Abundance values aligned with `taxonomies`.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by taxa just below the LCR with two columns:
        - 'ref_fracs'   : frequency-based fractions
        - 'abund_fracs': abundance-weighted fractions
    """
    if not lcr:
        return pd.DataFrame(columns=['ref_fracs', 'abund_fracs'], dtype=float)

    lcr_levels = lcr.split(';')
    lcr_depth = len(lcr_levels)

    # Extract taxa just below LCR
    next_level_taxa = []
    for tax in taxonomies['Taxonomy']:
        split_tax = tax.split(';')
        if len(split_tax) > lcr_depth:
            next_level_taxa.append(split_tax[lcr_depth])
        else:
            next_level_taxa.append('unclassified')

    # Count-based fractions
    ref_fracs = pd.Series(next_level_taxa).value_counts(normalize=True)

    # Abundance-weighted fractions
    df = pd.DataFrame({'taxon': next_level_taxa, 'abund': abunds})
    abund_sums = df.groupby('taxon')['abund'].sum()
    abund_fracs = abund_sums

    # Combine into one DataFrame
    result = pd.concat([ref_fracs, abund_fracs], axis=1)
    result.columns = ['ref_fracs', 'abund_fracs']

    return result


def build_gg_otu_id2tax_map(tree: Tree, tax_filepath: str, node: int, abunds: np.ndarray) -> dict:
    # interior_nodes = tree.find_interior_nodes()
    # if node not in interior_nodes:
    #     raise ValueError("Node must be an interior node!")

    table = load_tax_filepath(tax_filepath=tax_filepath)
    leaves = tree.find_leaves()
    otus = np.array([tree.get_name(i) for i in leaves]).astype(int)
    start = tree.find_subtree_start_indices()[node]
    size = tree.get_subtree_size()[node]

    taxonomies = table.loc[otus[start:start+size]]
    lcr = get_lowest_common_rank(taxonomies=taxonomies)
    fracs = get_relative_abundances_below_lcr(taxonomies, lcr, abunds[start:start+size])
    return {'LCR': lcr, 'Fracs': fracs}
