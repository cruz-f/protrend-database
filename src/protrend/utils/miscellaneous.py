import re
from collections import namedtuple
from datetime import datetime
from types import GeneratorType

from typing import Any

import numpy as np
import pandas as pd

from protrend.utils import Settings


WriteStack = namedtuple('WriteStack',
                        ['transformed', 'integrated', 'nodes', 'connected'],
                        defaults=[None, None, None, None])

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


def log_file_from_name(name: str) -> str:
    now = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    log_path = Settings.log_working_directory.joinpath(f'{name}_{now}.log')
    return log_path.as_posix()
