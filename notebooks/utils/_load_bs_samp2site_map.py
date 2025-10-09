# utils/_load_bs_samp2site_map.py
import pandas as pd


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