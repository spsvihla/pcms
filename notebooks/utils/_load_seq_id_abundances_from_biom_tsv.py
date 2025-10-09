# utils/_load_seq_id_abundances_from_biom_tsv.py
import pandas as pd


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
