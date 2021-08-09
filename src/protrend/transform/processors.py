import functools
import re
from collections import defaultdict
from typing import Callable, Any, List, Union

import pandas as pd
from w3lib.html import remove_tags

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
    for processor in processors:
        df[col] = df[col].map(processor)


def handle_nan(fn):
    """
    Decorator to handle NaN returns
    :param fn: Callable
    :return: the item or the results of the processor
    """
    @functools.wraps(fn)
    def wrapper(item):
        if item:
            return fn(item)
        return item

    return wrapper


@handle_nan
def operon_name(items: list) -> str:

    if items:

        names = defaultdict(int)

        for item in items:

            name = ''.join(letter for letter in item if letter.islower())

            if name:
                defaultdict[name] += 1

        name = items[0]
        best_score = 0

        for key, val in names.items():

            if val > best_score:
                best_score = val
                name = key

        return name

    return ''


@handle_nan
def genes_to_hash(items: list) -> str:
    items = sorted(items)
    return '_'.join(items)


@handle_nan
def str_join(items: list) -> str:
    return '_'.join(items)


@handle_nan
def null_to_nan(item: str) -> Union[None, str]:
    if item == 'null':
        return None
    return item


def nan_to_str(item: Any) -> str:
    if item:
        return item
    return ''


@handle_nan
def lstrip(item: str) -> str:
    return item.lstrip()


@handle_nan
def rstrip(item: str) -> str:
    return item.rstrip()


@handle_nan
def remove_white_space(item: str) -> str:
    return re.sub(white_space_pattern, repl='', string=item)


@handle_nan
def remove_multiple_white_space(item: str) -> str:
    return re.sub(multiple_white_space_pattern, repl=' ', string=item)


@handle_nan
def remove_regprecise_more(item: str) -> str:
    return re.sub(more_pattern, repl='', string=item)


@handle_nan
def remove_regprecise_more2(item: str) -> str:
    return re.sub(more2_pattern, repl='', string=item)


@handle_nan
def remove_pubmed(item: str) -> str:
    return re.sub(pubmed_pattern, repl='', string=item)


@handle_nan
def remove_html_tags(item: str) -> str:
    return remove_tags(item)


@handle_nan
def take_first(item: List[str]) -> str:
    if item:
        return item[0]

    return ''