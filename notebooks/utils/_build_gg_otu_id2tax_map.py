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


def get_relative_abundances_below_lcr(taxonomies: pd.DataFrame, lcr: str) -> pd.Series:
    """
    Given a set of taxonomy strings and the LCR, compute the relative abundances
    of the taxa at the rank just below the LCR.
    """
    if not lcr:
        return pd.Series(dtype=float)
    lcr_levels = lcr.split(';')
    lcr_depth = len(lcr_levels)

    next_level_taxa = []
    for tax in taxonomies['Taxonomy']:
        split_tax = tax.split(';')
        if len(split_tax) > lcr_depth:
            next_level_taxa.append(split_tax[lcr_depth])
        else:
            next_level_taxa.append('unclassified')

    counts = pd.Series(next_level_taxa).value_counts(normalize=True)
    return counts


def build_gg_otu_id2tax_map(tree: Tree, tax_filepath: str, node: int) -> dict:
    table = load_tax_filepath(tax_filepath=tax_filepath)
    leaves = tree.find_leaves(return_depths=False)
    otus = np.array([tree.get_name(i) for i in leaves]).astype(int)
    interior_nodes = tree.find_interior_nodes()
    subtree_sizes = tree.get_subtree_size()
    subtree_starts = tree.find_subtree_start_indices()

    if node not in interior_nodes:
        raise ValueError("Node must be an interior node!")
    clade = np.arange(subtree_starts[node], subtree_starts[node]+subtree_sizes[node])
    taxonomies = table.loc[otus[clade]]
    lcr = get_lowest_common_rank(taxonomies=taxonomies)
    rel_abund = get_relative_abundances_below_lcr(taxonomies, lcr)
    return {'LCR': lcr, 'RelAbundances': rel_abund}
