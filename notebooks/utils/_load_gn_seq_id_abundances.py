# utils/_load_gn_seq_id_abundances.py
import pandas as pd
from typing import Tuple


def load_gn_seq_id_abundances(
    table_filepath: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
