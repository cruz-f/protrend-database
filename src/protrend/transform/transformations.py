import re
from typing import Dict, Callable, List

import pandas as pd

from protrend.utils.processors import take_last


def create_input_value(df: pd.DataFrame, col: str) -> pd.DataFrame:
    return df.assign(input_value=df[col].copy())


def select_columns(df: pd.DataFrame, *columns: str) -> pd.DataFrame:
    df = df.copy()
    df = df.loc[:, list(columns)]
    return df


def merge_columns(df: pd.DataFrame, column: str, left: str, right: str) -> pd.DataFrame:
    df[column] = df[left].fillna(df[right])
    df = df.drop(columns=[left, right])
    return df


def merge_loci(df: pd.DataFrame, left_suffix: str, right_suffix: str) -> pd.DataFrame:
    # merge loci
    df = merge_columns(df=df, column='locus_tag',
                       left=f'locus_tag{left_suffix}', right=f'locus_tag{right_suffix}')
    df = df.dropna(subset=['locus_tag'])
    df = drop_empty_string(df, 'locus_tag')
    df = drop_duplicates(df=df, subset=['locus_tag'], perfect_match=True)
    return df


def concat_columns(df: pd.DataFrame, column: str, left: str, right: str) -> pd.DataFrame:
    df[column] = df[left] + df[right]
    df = df.drop(columns=[left, right])
    return df


def group_by(df: pd.DataFrame,
             column: str,
             aggregation: Dict[str, Callable],
             default: Callable = take_last) -> pd.DataFrame:
    agg = {col: default for col in df.columns if col != column}
    agg.update(aggregation)

    df = df.groupby(df[column]).aggregate(agg)
    df = df.reset_index()
    return df


def drop_empty_string(df: pd.DataFrame, *cols: str) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    df = df.copy()

    for col in cols:
        mask = df[col] != ''
        df = df[mask].copy()
        df = df.reset_index(drop=True)

    return df.reset_index(drop=True)


def drop_duplicates(df: pd.DataFrame,
                    subset: List[str],
                    perfect_match: bool = False,
                    preserve_nan: bool = True) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    if perfect_match and preserve_nan:

        df = df[(~df.duplicated(subset=subset)) | df[subset].isnull().any(axis=1)]

    elif perfect_match and not preserve_nan:

        df = df.drop_duplicates(subset=subset)

    elif not perfect_match and preserve_nan:

        for col in subset:
            df = df[(~df.duplicated(subset=[col])) | df[col].isnull()]

    elif not perfect_match and not preserve_nan:

        for col in subset:
            df = df.drop_duplicates(subset=[col])

    df = df.reset_index(drop=True)

    return df


def locus_tag_curation(df: pd.DataFrame) -> pd.DataFrame:
    bsu_pattern = r'^bsu[0-9]{5}$'
    ecoli_pattern = r'^b[0-9]{4}$'
    mtub_pattern = r'^rv[0-9]{4}$|^rv[0-9]{4}c$'
    cg_pattern = r'^cg[0-9]{4}$'

    for _, row in df.iterrows():
        tax = row['taxonomy']
        locus_tag = row['locus_tag']
        locus_tag_to_set = None
        synonyms = row.get('synonyms', [])

        if tax == '224308' and not re.match(bsu_pattern, locus_tag, re.IGNORECASE):

            for synonym in synonyms:
                if re.match(bsu_pattern, synonym, re.IGNORECASE):
                    locus_tag_to_set = synonym
                    break

        elif tax == '511145' and not re.match(ecoli_pattern, locus_tag, re.IGNORECASE):
            for synonym in synonyms:
                if re.match(ecoli_pattern, synonym, re.IGNORECASE):
                    locus_tag_to_set = synonym
                    break

        elif tax == '83332' and not re.match(mtub_pattern, locus_tag, re.IGNORECASE):
            for synonym in synonyms:
                if re.match(mtub_pattern, synonym, re.IGNORECASE):
                    locus_tag_to_set = synonym
                    break

        elif tax == '196627' and not re.match(cg_pattern, locus_tag, re.IGNORECASE):
            for synonym in synonyms:
                if re.match(cg_pattern, synonym, re.IGNORECASE):
                    locus_tag_to_set = synonym
                    break

        else:
            continue

        row['locus_tag'] = locus_tag_to_set

    df = df.dropna(subset=['locus_tag'])
    df = drop_empty_string(df, 'locus_tag')
    return df
