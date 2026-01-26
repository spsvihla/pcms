from typing import Optional
from numpy.typing import NDArray
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


def get_lowest_common_rank(taxa: pd.DataFrame) -> str:
    tax_split = [row['Taxonomy'].split(';') for _, row in taxa.iterrows()]
    ranks_by_level = list(zip(*tax_split))
    lca = []
    for level in ranks_by_level:
        if len(set(level)) == 1 and level[0]:
            lca.append(level[0])
        else:
            break
    return ';'.join(lca)


def get_relative_abundances_below_lcr(
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
    table = load_tax_filepath(tax_filepath=tax_filepath)
    leaves = tree.find_leaves()
    otus = np.array([tree.get_name(i) for i in leaves]).astype(int)
    start = tree.find_subtree_start_indices()[node]
    size = tree.get_subtree_size()[node]

    taxa = table.loc[otus[start:start+size]]
    lcr = get_lowest_common_rank(taxa=taxa)
    if b_abunds is None:
        fracs = get_relative_abundances_below_lcr(taxa, lcr, a_abunds[start:start+size])
    else:
        fracs = get_relative_abundances_below_lcr(taxa, lcr, a_abunds[start:start+size], b_abunds[start:start+size])
    
    return {'lcr': lcr, 'fracs': fracs}
