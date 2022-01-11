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


def merge_loci(self, df: pd.DataFrame, left_suffix: str, right_suffix: str) -> pd.DataFrame:
    # merge loci
    df = self.merge_columns(df=df, column='locus_tag',
                            left=f'locus_tag{left_suffix}', right=f'locus_tag{right_suffix}')
    df = df.dropna(subset=['locus_tag'])
    df = self.drop_empty_string(df, 'locus_tag')
    df = self.drop_duplicates(df=df, subset=['locus_tag'], perfect_match=True)
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
    n_rows, _ = df.shape

    mask = pd.Series([False] * n_rows)

    cols_mask = [df[col] != '' for col in cols]

    if cols_mask:
        mask = pd.concat(cols_mask, axis=1)
        mask = mask.all(axis=1)

    new_df = df[mask].copy()
    return new_df


def drop_duplicates(df: pd.DataFrame,
                    subset: List[str],
                    perfect_match: bool = False,
                    preserve_nan: bool = True) -> pd.DataFrame:
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
