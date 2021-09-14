from typing import Set, Dict, Callable

import pandas as pd


def read_from_stack(stack: Dict[str, str],
                    file: str,
                    default_columns: Set[str],
                    reader: Callable,
                    **kwargs) -> pd.DataFrame:

    try:
        file_path = stack[file]
        return reader(file_path, **kwargs)

    except (OSError, KeyError):

        default_columns = list(default_columns)
        return pd.DataFrame(columns=default_columns)
