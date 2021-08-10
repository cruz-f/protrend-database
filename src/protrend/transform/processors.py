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
def remove_ellipsis(item: str) -> str:
    if item.endswith('...'):
        return item[:-3]

    return item


@handle_nan
def upper_case(item: str) -> str:
    return item.upper()


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


def operon_strand(previous_strand: str = None,
                  current_strand: str = None,
                  default_strand: str = 'NA') -> str:
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


def operon_left_position(strand: str,
                         previous_left: int = None,
                         current_left: int = None,
                         default_left: Union[str, int] = 'NA') -> Union[str, int]:
    left = default_left

    if previous_left is None:
        previous_left = 0

    if current_left is None:
        current_left = 0

    if strand in ('forward', 'NA'):

        if current_left < previous_left:
            left = current_left

        elif previous_left < current_left:
            left = previous_left

        elif previous_left == current_left and current_left != 0:
            left = previous_left

    elif strand == 'reverse':

        if current_left > previous_left:
            left = current_left

        elif previous_left > current_left:
            left = previous_left

        elif previous_left == current_left and current_left != 0:
            left = previous_left

    return left


def operon_right_position(strand: str,
                          previous_right: int = None,
                          current_right: int = None,
                          default_right: Union[str, int] = 'NA') -> Union[str, int]:
    right = default_right

    if previous_right is None:
        previous_right = 0

    if current_right is None:
        current_right = 0

    if strand in ('forward', 'NA'):

        if current_right > previous_right:
            right = current_right

        elif previous_right > current_right:
            right = previous_right

        elif previous_right == current_right and current_right != 0:
            right = previous_right

    elif strand == 'reverse':

        if current_right < previous_right:
            right = current_right

        elif previous_right < current_right:
            right = previous_right

        elif previous_right == current_right and current_right != 0:
            right = previous_right

    return right


def tfbs_left_position(strand: str,
                       gene_position: Union[str, int, None],
                       gene_relative_position: int,
                       default: Union[str, int] = 'NA'):
    if gene_position == 'NA' or gene_position is None:
        return default

    if strand in ('forward', 'NA'):
        return gene_position + gene_relative_position

    elif strand == 'reverse':
        return gene_position - gene_relative_position

    return default


def tfbs_right_position(strand: str,
                        gene_position: Union[str, int],
                        gene_relative_position: int,
                        tfbs_length: int,
                        default: Union[str, int] = 'NA'):

    if gene_position == 'NA' or gene_position is None:
        return default

    if strand in ('forward', 'NA'):
        return gene_position + (gene_relative_position + tfbs_length)

    elif strand == 'reverse':
        return gene_position - (gene_relative_position + tfbs_length)

    return default
