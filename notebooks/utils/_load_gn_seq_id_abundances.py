# utils/_load_gn_seq_id_abundances.py
import pandas as pd
from typing import Tuple


def load_gn_seq_id_abundances(
    data_filepath: str,
    metadata_filepath: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load sequence abundance data and associated metadata from TSV files.

    Parameters
    ----------
    data_filepath : str
        Path to the TSV file containing sequence abundance data.
    metadata_filepath : str
        Path to the TSV file containing sample metadata.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        A tuple containing:
        - data: DataFrame of abundance values (rows are sequence IDs, columns are samples)
        - metadata: DataFrame of metadata (indexed by sample or with sample identifiers as a column)
    """
    try:
        data = pd.read_csv(data_filepath, sep='\t', engine='c')
        metadata = pd.read_csv(metadata_filepath, sep='\t', engine='c')
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Missing file: {e.filename}") from e
    except pd.errors.ParserError as e:
        raise RuntimeError(f"Parsing error: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Unexpected error while loading files: {e}") from e

    return data, metadata
