import os
from pathlib import Path
from typing import Union, IO, AnyStr

import pandas as pd


def read_csv(filepath_or_buffer: Union[str, Path, IO[AnyStr]], **kwargs) -> pd.DataFrame:

    if not os.path.exists(filepath_or_buffer):
        raise OSError(f'Invalid file path {filepath_or_buffer}')

    return pd.read_csv(filepath_or_buffer, **kwargs)