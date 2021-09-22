import re
from collections import defaultdict
from typing import Callable, Any, List, Union, Sequence, Set

import pandas as pd
from w3lib.html import remove_tags

from protrend.utils.miscellaneous import is_null

more_pattern = re.compile(r'\s\smore\s\s\s')
more2_pattern = re.compile(r'\s\smore\s')
white_space_pattern = re.compile(r'\s')
multiple_white_space_pattern = re.compile(r' +')
pubmed_pattern = re.compile(r'\[[0-9]*]|\[[0-9]*,|\s[0-9]*,|\s[0-9]*]')
pubmed_pattern2 = re.compile(r'\[pmid::[0-9]*]|\[pmid: [0-9]*]|\[pmid:: [0-9]*]|\[pmid:[0-9]*]|\[ pmid: : [0-9]* ]|\[pmid: : [0-9]*]|\[pmid::[0-9]* ]|\[ pmid::[0-9]*]|\[pmid::[0-9]* .]')


def apply_processors(df: pd.DataFrame, **processors: Union[Callable, List[Callable]]) -> pd.DataFrame:
    """
    Helping function to apply processors over a pandas DataFrame column

    :param processors: multiple processors to be applied over DataFrame columns
    :param df: DataFrame that will be changed by multiple processors applied over a single column

    :return: it returns a new DataFrame object

    """

    handle_nan_processors = (null_to_str, null_to_none, to_nan, to_list_nan)

    new_columns = {}

    for column, processor in processors.items():

        processor = to_list(processor)
        ds = df[column].copy()

        for fn in processor:

            if fn in handle_nan_processors:
                ds = ds.map(fn)

            else:
                ds = ds.map(fn, na_action='ignore')

        new_columns[column] = ds

    df = df.assign(**new_columns)

    return df


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

    item_to_int = to_int(item)

    try:
        return str(item_to_int)
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


def remove_pubmed2(item: str) -> str:
    return re.sub(pubmed_pattern2, repl='', string=item)


def remove_html_tags(item: str) -> str:
    return remove_tags(item)


def parse_collectf_pubmed(item: List[str]) -> List[str]:

    escapes = ''.join([chr(char) for char in range(1, 32)])
    translator = str.maketrans('', '', escapes)

    all_pmids = set()
    for string in item:

        to_remove = string.replace(' ', '').replace('.', '')

        if to_remove == '':
            continue

        if string.startswith('>a') or string.startswith('>c'):
            continue

        string = string.translate(translator)

        pmids1 = re.findall(pubmed_pattern, string)
        pmids2 = re.findall(pubmed_pattern2, string)
        pmids = pmids1 + pmids2
        for pmid in pmids:
            pmid = ''.join(l for l in pmid if l.isdigit())
            if pmid:
                all_pmids.add(pmid)

    return list(all_pmids)


def parse_collectf_description(item: List[str]) -> str:

    escapes = ''.join([chr(char) for char in range(1, 32)])
    translator = str.maketrans('', '', escapes)

    to_concat = []
    for string in item:

        to_remove = string.replace(' ', '').replace('.', '')

        if to_remove == '':
            continue

        if string.startswith('>a') or string.startswith('>c'):
            continue

        string = string.translate(translator)
        string = remove_pubmed2(string)
        string = remove_pubmed(string)
        to_concat.append(string)

    return ''.join(to_concat)


def regulatory_effect(items: Set[str]) -> Union[None, str]:
    new_items = {item.replace(' ', '') for item in items if not is_null(item)}

    if len(new_items) == 0:
        return

    if len(new_items) > 1:
        return 'dual'

    if 'repressor' in new_items:
        return 'repression'

    if 'activator' in new_items:
        return 'activation'

    return 'dual'


def regulatory_effect_collectf(item: str) -> Union[None, str]:

    if is_null(item):
        return

    if item.lower() == 'rep':
        return 'repression'

    if item.lower() == 'act':
        return 'activation'

    return


def parse_effector_name_regulondb(item: str) -> str:
    if is_null(item):
        return

    return item.replace('&', '').replace(';', '')
