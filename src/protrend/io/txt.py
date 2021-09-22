import os
from pathlib import Path
from typing import Union, IO, AnyStr

import pandas as pd


def read_txt(filepath_or_buffer: Union[str, Path, IO[AnyStr]], **kwargs) -> pd.DataFrame:
    """
    Wrapper for pandas TXT reader.
    This method is optimized to read TXT files exported by the transformers,
    and is mostly used in the transformers and loaders.


    :param filepath_or_buffer: file path as string or file handler
    :param kwargs: additional kwargs for pd.read_csv API

    :return: pandas DataFrame out of the TXT file
    """

    if not os.path.exists(filepath_or_buffer):
        raise OSError(f'Invalid file path {filepath_or_buffer}')

    if 'sep' not in kwargs:
        kwargs['sep'] = '\t'

    return pd.read_csv(filepath_or_buffer, **kwargs)
