from typing import Dict, Callable, Sequence

import pandas as pd

from protrend.utils import MultiStack


def _default_df(columns) -> pd.DataFrame:
    columns = list(columns)
    return pd.DataFrame(columns=columns)


def read_from_stack(stack: Dict[str, str],
                    key: str,
                    columns: Sequence[str],
                    reader: Callable) -> pd.DataFrame:

    try:
        file_path = stack[key]
        return reader(file_path)

    except (OSError, KeyError):

        return _default_df(columns)


def read_from_multi_stack(stack: Dict[str, MultiStack],
                          key: str,
                          columns: Sequence[str]) -> pd.DataFrame:

    try:
        multi_stack = stack[key]

        dfs = []
        for file, taxon, source, reader in zip(multi_stack.stack, multi_stack.taxa,
                                               multi_stack.source, multi_stack.reader):
            df = reader(file)
            df = df.assign(taxonomy=taxon, source=source)
            dfs.append(df)

        return pd.concat(dfs)

    except (OSError, KeyError):

        df = _default_df(columns)

        default_cols = {}
        if 'taxonomy' not in df.columns:
            default_cols['taxonomy'] = None
        if 'source' not in df.columns:
            default_cols['source'] = None

        df = df.assign(**default_cols)
        return df
