from typing import Dict, Callable, Sequence

import pandas as pd

from .json import read_json_frame
from protrend.utils import build_file_path


def read(source: str,
         version: str,
         file: str,
         reader: Callable,
         kwargs: Dict = None,
         default: pd.DataFrame = None) -> pd.DataFrame:
    if not kwargs:
        kwargs = {}

    file_path = build_file_path(source=source, version=version, file=file)

    try:
        return reader(file_path, **kwargs)

    except OSError:

        return default


def read_source(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='integrated_source.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))


def read_organism(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='integrated_organism.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))


def read_regulator(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='integrated_regulator.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))


def read_gene(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='integrated_gene.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))


def read_tfbs(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='integrated_tfbs.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))


def read_effector(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='integrated_effector.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))


def read_rfam(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='integrated_regulatoryfamily.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))


def read_operon(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='integrated_operon.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))


def read_promoters(source: str, version: str, columns: Sequence[str]) -> pd.DataFrame:
    return read(source=source, version=version, file='transformed_promoters.json',
                reader=read_json_frame, default=pd.DataFrame(columns=columns))
