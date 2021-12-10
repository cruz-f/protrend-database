import os
import re
from collections import namedtuple
from types import GeneratorType
from typing import Any, Dict, List

import numpy as np
import pandas as pd

from protrend.utils import Settings

WriteStack = namedtuple('WriteStack',
                        ['transformed', 'integrated', 'nodes', 'connected'],
                        defaults=[None, None, None, None])


def build_stack(source: str,
                version: str,
                stack_to_load: Dict[str, str],
                sa: bool = True,
                dl: bool = True) -> Dict[str, str]:
    loaded_stack = {}

    for key, file in stack_to_load.items():

        if sa:
            file = os.path.join(Settings.staging_area, source, version, file)
            loaded_stack[key] = file

        if dl:
            file = os.path.join(Settings.data_lake, source, version, file)
            loaded_stack[key] = file

    return loaded_stack


def build_load_stack(source: str,
                     version: str,
                     stack_to_load: List[str]) -> List[str]:
    loaded_stack = []

    for file in stack_to_load:
        file = os.path.join(Settings.data_lake, source, version, f'{file}.json')
        loaded_stack.append(file)

    return loaded_stack


# CamelCase to snake_case
camel_case_pattern = re.compile(r'(?<!^)(?=[A-Z])')


def convert_to_snake_case(item: str):
    return camel_case_pattern.sub('_', item).lower()


def args_length(*args):
    size = 0

    for param in args:

        if param is not None:
            param_size = len(param)
            if param_size > size:
                size = param_size

    return size


def scale_arg(arg, size):
    if arg is None:
        return [None] * size

    elif len(arg) == 0:
        return [None] * size

    elif len(arg) == 1:
        return arg * size

    elif len(arg) != size:
        raise ValueError(f'Invalid input size of {len(arg)}')

    return arg


def is_null(obj: Any) -> bool:
    # booleans
    if isinstance(obj, bool):
        return not obj

    # integers, floats, etc
    if isinstance(obj, (int, float)):
        return pd.isnull(obj)

    # numpy arrays
    if isinstance(obj, np.ndarray):
        return not obj.any()

    # pandas series and frames
    if isinstance(obj, (pd.DataFrame, pd.Series)):
        return obj.empty

    # python built-ins
    if isinstance(obj, (range, list, tuple, set, dict, frozenset, str)):
        return len(obj) == 0

    # python generator built-ins
    if isinstance(obj, GeneratorType):
        return False

    # pandas check for nan or null
    if pd.isnull(obj):
        return True

    return False
