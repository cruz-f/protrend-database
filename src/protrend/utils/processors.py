import re
from statistics import mode, StatisticsError
from typing import Callable, Any, List, Union, Sequence

import numpy as np
import pandas as pd
from w3lib.html import remove_tags

from protrend.utils import SetList
from protrend.utils.constants import UNKNOWN, ACTIVATION, REPRESSION, DUAL
from protrend.utils.miscellaneous import is_null

more_pattern = re.compile(r'\s\smore\s\s\s')
more2_pattern = re.compile(r'\s\smore\s')
white_space_pattern = re.compile(r'\s')
multiple_white_space_pattern = re.compile(r' +')
pubmed_pattern = re.compile(r'\[[0-9]*]|\[[0-9]*,|\s[0-9]*,|\s[0-9]*]')
pubmed_pattern2 = re.compile(r'\[pmid::[0-9]*]|\[pmid: [0-9]*]|\[pmid:: [0-9]*]|\[pmid:[0-9]*]|\[ pmid: : [0-9]* ]|\[pmid: : [0-9]*]|\[pmid::[0-9]* ]|\[ pmid::[0-9]*]|\[pmid::[0-9]* .]')
pubmed_pattern3 = re.compile(r'\[PMID::[0-9]*]|\[PMID: [0-9]*]|\[PMID:: [0-9]*]|\[PMID:[0-9]*]|\[ PMID: : [0-9]* ]|\[PMID: : [0-9]*]|\[PMID::[0-9]* ]|\[ PMID::[0-9]*]|\[PMID::[0-9]* .]')


def apply_processors(df: pd.DataFrame, **processors: Union[Callable, List[Callable]]) -> pd.DataFrame:
    """
    Helping function to apply processors over a pandas DataFrame column

    :param processors: multiple processors to be applied over DataFrame columns
    :param df: DataFrame that will be changed by multiple processors applied over a single column

    :return: it returns a new DataFrame object

    """

    handle_nan_processors = (null_to_str, null_to_none, to_nan, to_list_nan, flatten_set_list_nan, protrend_hash)
    regulatory_effect_processors = (regulatory_effect_abasy, regulatory_effect_collectf, regulatory_effect_coryneregnet,
                                    regulatory_effect_dbtbs, regulatory_effect_literature, regulatory_effect_literature,
                                    regulatory_effect_regprecise, regulatory_effect_regulondb)

    new_columns = {}

    for column, processor in processors.items():

        processor = to_list(processor)
        ds = df[column].copy()

        for fn in processor:

            if fn in handle_nan_processors or fn in regulatory_effect_processors:
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


def to_int(item: Any) -> int:
    if isinstance(item, int):
        return item

    if is_null(item):
        return item

    try:
        return int(item)
    except (ValueError, TypeError):
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


def to_list(item: Any) -> list:
    if isinstance(item, str):
        return [item]

    try:
        iterator = iter(item)
    except TypeError:
        iterator = iter([item])

    return list(iterator)


def to_set_list(item: Any) -> SetList:
    if isinstance(item, str):
        return SetList([item])

    try:
        iterator = iter(item)
    except TypeError:
        iterator = iter([item])

    return SetList(iterator)


def to_nan(item: Any) -> Union[None, Any]:
    if is_null(item):
        return

    return item


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


def flatten_list(items):
    try:
        res = []
        for item in items:
            try:
                for sub_item in item:
                    res.append(sub_item)

            except TypeError:
                continue

        return res

    except TypeError:
        return []


def flatten_set_list_nan(items):
    flatten_lst = flatten_list(items=items)
    return SetList(flatten_lst)


def take_last(items: Sequence[Any]) -> Any:
    if isinstance(items, pd.Series):
        try:
            return items.iloc[-1]

        except IndexError:
            return

    try:
        return items[-1]

    except IndexError:
        return


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


def split_semi_colon(item: str) -> List[str]:
    items = item.split(sep=';')
    return [i.rstrip().lstrip() for i in items]


def upper_case(item: str) -> str:
    return item.upper()


def lower_case(item: str) -> str:
    return item.lower()


def strand_mode(item):
    if is_null(item):
        return UNKNOWN

    try:
        m = mode(item)

        if is_null(m):
            return UNKNOWN

        return m

    except StatisticsError:

        for sub_item in item:
            return sub_item


def start_forward(item):
    if is_null(item):
        return

    item = to_list(item)

    x = np.array(item, dtype=np.float64)
    return np.nanmin(x)


def start_reverse(item):
    if is_null(item):
        return

    item = to_list(item)

    x = np.array(item, dtype=np.float64)
    return np.nanmax(x)


def protrend_hash(items: List[str]) -> Union[None, str]:
    if is_null(items):
        return

    items = [to_str(item) for item in items]

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


def remove_pubmed3(item: str) -> str:
    return re.sub(pubmed_pattern3, repl='', string=item)


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
            pmid = ''.join(digit for digit in pmid if digit.isdigit())
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
        string = remove_pubmed3(string)
        string = remove_pubmed(string)
        to_concat.append(string)

    return ''.join(to_concat)


def regulatory_effect_regprecise(item: str) -> str:
    if is_null(item):
        return UNKNOWN

    item = item.lstrip().rstrip()

    if 'repressor' == item:
        return REPRESSION

    if 'activator' == item:
        return ACTIVATION

    if 'dual' == item:
        return DUAL

    if 'activator' in item and 'repressor' in item:
        return DUAL

    return UNKNOWN


def regulatory_effect_collectf(item: str) -> str:
    if is_null(item):
        return UNKNOWN

    if item.lower() == 'rep':
        return REPRESSION

    if item.lower() == 'act':
        return ACTIVATION

    return UNKNOWN


def regulatory_effect_regulondb(item: str) -> str:
    if is_null(item):
        return UNKNOWN

    if item.lower() == '-' or item.lower() == 'repressor':
        return REPRESSION

    if item.lower() == '+' or item.lower() == 'activator':
        return ACTIVATION

    if item.lower() == '+-' or item.lower() == DUAL:
        return DUAL

    return UNKNOWN


def regulatory_effect_dbtbs(item: str) -> str:
    if is_null(item):
        return UNKNOWN

    item = item.rstrip().lstrip()

    if item.lower() == 'negative':
        return REPRESSION

    if item.lower() == 'positive':
        return ACTIVATION

    if item.lower() == 'promoter' or item.lower() == 'promorter':
        return ACTIVATION

    if item.lower() == 'pos/neg':
        return DUAL

    if item.lower() == 'positive (siga promoter); negative (sige promoter)':
        return DUAL

    if item.lower() == 'negative (sige promoter)':
        return REPRESSION

    if item.lower() == 'negative at high mn(ii) concentration':
        return REPRESSION

    if item.lower() == 'positive at low mn(ii) concentration':
        return ACTIVATION

    return UNKNOWN


def regulatory_effect_coryneregnet(item: str) -> str:
    if is_null(item):
        return UNKNOWN

    if item.lower() == 'repressor':
        return REPRESSION

    if item.lower() == 'activator':
        return ACTIVATION

    return UNKNOWN


def regulatory_effect_abasy(item: str) -> str:
    if is_null(item):
        return UNKNOWN

    if item.lower() == '-':
        return REPRESSION

    if item.lower() == '+':
        return ACTIVATION

    if item.lower() == '?':
        return UNKNOWN

    return UNKNOWN


def _parse_regulatory_effect_literature_bsub(item: str):
    signs = item.split('|')

    is_positive = False
    is_negative = False
    for sign in signs:

        if sign.rstrip().lstrip().lower() == '-1':
            is_negative = True

        if sign.rstrip().lstrip().lower() == '1':
            is_positive = True

    if is_positive and is_negative:
        return DUAL

    if is_positive:
        return ACTIVATION

    if is_negative:
        return REPRESSION

    return UNKNOWN


def regulatory_effect_literature(item: str) -> str:
    if is_null(item):
        return UNKNOWN

    item = item.rstrip().lstrip().lower()

    if '|' in item:
        return _parse_regulatory_effect_literature_bsub(item)

    if item == '+' or item == 'ind':
        return ACTIVATION

    if item == '-' or item == 'rep':
        return REPRESSION

    if item == '?' or item == '.':
        return UNKNOWN

    if item == 'd':
        return DUAL

    return UNKNOWN


def parse_effector_name_regulondb(item: str) -> Union[None, str]:
    if is_null(item):
        return

    return item.replace('&', '').replace(';', '')
