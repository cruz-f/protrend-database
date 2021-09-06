import re
from collections import defaultdict
from typing import Callable, Any, List, Union, Sequence

import numpy as np
import pandas as pd
from w3lib.html import remove_tags

from protrend.utils.miscellaneous import is_null

more_pattern = re.compile(r'\s\smore\s\s\s')
more2_pattern = re.compile(r'\s\smore\s')
white_space_pattern = re.compile(r'\s')
multiple_white_space_pattern = re.compile(r' +')
pubmed_pattern = re.compile(r'\[[0-9]*]|\[[0-9]*,|\s[0-9]*,|\s[0-9]*]')


def apply_processors(*processors: Callable, df: pd.DataFrame, col: str) -> None:
    """
    Helping function to apply processors over a pandas DataFrame column

    :param processors: multiple processors to be applied over a DataFrame column
    :param df: DataFrame that will be changed by multiple processors applied over a single column
    :param col: DataFrame column

    :return: Nothing to return, it just changes the DataFrame object

    """

    handle_nan_processors = (null_to_str, null_to_none, to_nan, to_list_nan)

    for processor in processors:
        if processor in handle_nan_processors:
            df[col] = df[col].map(processor)

        else:
            df[col] = df[col].map(processor, na_action='ignore')


def to_str(item: Any) -> str:
    if isinstance(item, str):
        return item

    if is_null(item):
        return item

    try:
        return str(item)
    except (ValueError, TypeError):
        return item


def to_nan(item: Any) -> Union[None, Any]:
    if is_null(item):
        return None

    return item


def to_int_str(item: Any) -> str:
    if isinstance(item, str):
        return item

    if is_null(item):
        return item

    item = to_int(item)

    try:
        return str(item)
    except (ValueError, TypeError):
        return item


def to_int(item: Any) -> int:
    if isinstance(item, int):
        return item

    if is_null(item):
        return item

    try:
        return int(item)
    except (ValueError, TypeError):
        return item


def to_set(item: Any) -> set:
    if isinstance(item, str):
        return {item}

    try:
        iterator = iter(item)
    except TypeError:
        iterator = iter([item])

    return set(iterator)


def to_list(item: Any) -> list:
    if isinstance(item, str):
        return [item]

    try:
        iterator = iter(item)
    except TypeError:
        iterator = iter([item])

    return list(iterator)


def to_list_nan(item: Any) -> list:
    if isinstance(item, str):
        return [item]

    if is_null(item):
        return []

    try:
        iterator = iter(item)
    except TypeError:
        iterator = iter([item])

    return list(iterator)


def flatten_set(items):
    return {i for arg in items for i in arg}


def flatten_list(items):
    return [i for arg in items for i in arg]


def take_last(items: Sequence[Any]) -> Any:

    if isinstance(items, pd.Series):
        try:
            return items.iloc[-1]

        except IndexError:
            return None

    try:
        return items[-1]

    except IndexError:
        return None


def take_first(items: Sequence[Any]) -> Any:

    if isinstance(items, pd.Series):
        try:
            return items.iloc[0]

        except IndexError:
            return None

    try:
        return items[0]

    except IndexError:
        return None


def null_to_none(item: Any) -> Union[Any, None]:

    if is_null(item):
        return None

    return item


def null_to_str(item: Any) -> str:
    if is_null(item):
        return ''
    return item


def split_str(item: str) -> List[str]:
    return item.split(sep=' ')


def upper_case(item: str) -> str:
    return item.upper()


def lower_case(item: str) -> str:
    return item.lower()


def operon_name(items: List[str]) -> str:
    if items:

        names = defaultdict(int)

        for item in items:

            name = ''.join(letter for letter in item if letter.islower())

            if len(name) > 2:
                names[name] += 1

        name = ''
        best_score = 0

        for key, val in names.items():

            if val > best_score:
                best_score = val
                name = key

        return name

    return ''


def genes_to_hash(items: List[str]) -> str:
    items = sorted(items)
    return '_'.join(items)


def str_join(items: List[str]) -> str:
    return '_'.join(items)


def lstrip(item: str) -> str:
    return item.lstrip()


def rstrip(item: str) -> str:
    return item.rstrip()


def remove_ellipsis(item: str) -> str:
    if item.endswith('...'):
        return item[:-3]

    return item


def remove_white_space(item: str) -> str:
    return re.sub(white_space_pattern, repl='', string=item)


def remove_multiple_white_space(item: str) -> str:
    return re.sub(multiple_white_space_pattern, repl=' ', string=item)


def remove_regprecise_more(item: str) -> str:
    return re.sub(more_pattern, repl='', string=item)


def remove_regprecise_more2(item: str) -> str:
    return re.sub(more2_pattern, repl='', string=item)


def remove_pubmed(item: str) -> str:
    return re.sub(pubmed_pattern, repl='', string=item)


def remove_html_tags(item: str) -> str:
    return remove_tags(item)


def operon_strand(previous_strand: str = None,
                  current_strand: str = None,
                  default_strand: str = None) -> Union[None, str]:
    strand = None

    if previous_strand is None:
        previous_strand = ''

    if current_strand is None:
        current_strand = ''

    if current_strand in ('forward', 'reverse'):
        strand = current_strand

    elif previous_strand in ('forward', 'reverse'):
        strand = previous_strand

    if strand is None:
        strand = default_strand

    return strand


def operon_left_position(strand: Union[None, str],
                         previous_left: int = None,
                         current_left: int = None,
                         default_left: Union[None, int] = None) -> Union[None, int]:
    left = default_left

    if strand is None:
        strand = 'forward'

    if strand == 'forward':

        if previous_left is None:
            previous_left = np.inf

        if current_left is None:
            current_left = np.inf

        if current_left < previous_left:
            left = current_left

        elif previous_left < current_left:
            left = previous_left

        elif previous_left == current_left and current_left != np.inf:
            left = current_left

    elif strand == 'reverse':

        if previous_left is None:
            previous_left = -np.inf

        if current_left is None:
            current_left = -np.inf

        if current_left > previous_left:
            left = current_left

        elif previous_left > current_left:
            left = previous_left

        elif previous_left == current_left and current_left != 0:
            left = current_left

    return left


def operon_right_position(strand: str,
                          previous_right: int = None,
                          current_right: int = None,
                          default_right: Union[None, int] = None) -> Union[None, int]:
    right = default_right

    if strand is None:
        strand = 'forward'

    if strand == 'forward':

        if previous_right is None:
            previous_right = np.inf

        if current_right is None:
            current_right = np.inf

        if current_right > previous_right:
            right = current_right

        elif previous_right > current_right:
            right = previous_right

        elif previous_right == current_right and current_right != 0:
            right = current_right

    elif strand == 'reverse':

        if previous_right is None:
            previous_right = -np.inf

        if current_right is None:
            current_right = -np.inf

        if current_right < previous_right:
            right = current_right

        elif previous_right < current_right:
            right = previous_right

        elif previous_right == current_right and current_right != 0:
            right = previous_right

    return right


def tfbs_left_position(strand: str,
                       gene_position: Union[int, None],
                       gene_relative_position: int,
                       default: Union[None, int] = None) -> Union[None, int]:

    if strand is None:
        strand = 'forward'

    if gene_position is None:
        return default

    if strand == 'forward':
        return gene_position + gene_relative_position

    elif strand == 'reverse':
        return gene_position - gene_relative_position

    return default


def tfbs_right_position(strand: str,
                        gene_position: Union[int, None],
                        gene_relative_position: int,
                        tfbs_length: int,
                        default: Union[None, int] = None) -> Union[None, int]:

    if strand is None:
        strand = 'forward'

    if gene_position is None:
        return default

    if strand == 'forward':
        return gene_position + (gene_relative_position + tfbs_length)

    elif strand == 'reverse':
        return gene_position - (gene_relative_position + tfbs_length)

    return default
