import json
import os
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from protrend.utils import SetList


def read_json(file_path: str) -> dict:
    if not os.path.exists(file_path):
        raise OSError(f'Invalid file path {file_path}')

    with open(file_path, 'r') as handle:
        return json.load(handle)


def write_json(file_path: str, content: dict):
    with open(file_path, 'w') as handle:
        return json.dump(content, handle)


def read_json_lines(file_path: str) -> pd.DataFrame:
    if not os.path.exists(file_path):
        raise OSError(f'Invalid file path {file_path}')

    lines = []

    with open(file_path, 'r') as file:

        for line in file:
            line = line.rstrip()

            if line:
                line_dict = json.loads(line)
                lines.append(line_dict)

    return pd.DataFrame(lines)


def read_json_frame(file_path: Union[Path, str], **kwargs) -> pd.DataFrame:
    """
        Wrapper for pandas JSON reader.
        This method is optimized to read JSON files exported by the transformers,
        and is mostly used in the transformers and loaders.


        :param file_path: file path as string
        :param kwargs: additional kwargs for pd.read_json API

        :return: pandas DataFrame out of the JSON file
    """

    if not os.path.exists(file_path):
        raise OSError(f'Invalid file path {file_path}')

    return pd.read_json(file_path, **kwargs)


def _set_list_serializer(o: SetList) -> list:
    return o.data


def write_json_frame(file_path: str, df: pd.DataFrame, **kwargs) -> Optional[str]:

    if 'default_handler' not in kwargs:
        kwargs['default_handler'] = _set_list_serializer

    return df.to_json(file_path, **kwargs)
