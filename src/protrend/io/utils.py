from typing import Set, Union, TYPE_CHECKING

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines


if TYPE_CHECKING:
    from protrend.transform.connector import Connector
    from protrend.transform.transformer import Transformer


def read_from_stack(tl: Union['Transformer', 'Connector'],
                    file: str,
                    json: bool,
                    default_columns: Set[str]) -> pd.DataFrame:

    file_path = None

    if hasattr(tl, 'transform_stack'):

        file_path = tl.transform_stack.get(file)

    elif hasattr(tl, 'connect_stack'):

        file_path = tl.connect_stack.get(file)

    if file_path:
        if json:
            df = read_json_lines(file_path)

        else:
            df = read_csv(file_path)

    else:
        default_columns = list(default_columns)
        df = pd.DataFrame(columns=default_columns)

    return df