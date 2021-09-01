from typing import Set, Dict, Callable

import pandas as pd


def read_from_stack(stack: Dict[str, str],
                    file: str,
                    default_columns: Set[str],
                    reader: Callable,
                    **kwargs) -> pd.DataFrame:

    if file in stack:

        file_path = stack[file]

        df = reader(file_path, **kwargs)

    else:
        default_columns = list(default_columns)
        df = pd.DataFrame(columns=default_columns)

    return df
