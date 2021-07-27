import os
from abc import ABCMeta, abstractmethod
from typing import Dict, Tuple

import pandas as pd

from protrend.io.json import read_json_lines
from protrend.model.node import Node
from protrend.utils.settings import STAGING_AREA_PATH, DATA_LAKE_PATH


class Transformer(metaclass=ABCMeta):
    """
    Transformer interface.

    The following methods must be implemented to set up a transformer for each node
    and relationship

    A transformer is responsible for reading, processing, integrating and writing node and relationships files.

    A transformer starts with data from the staging area and ends with structured nodes and relationships.
    """
    node: Node = None

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 **files: Dict[str, str]):

        """
        The transform object contains and uses transformation procedures for a given neomodel node entity.
        This object is responsible for reading several files and digest/process them into node instances.

        A given source/database contained in the staging area must be provided
        together with the files required for the transformation to take place.

        An attribute will be created for each key-value pair in the files' mapping-like object.

        A pandas DataFrame is the main engine to read, load, process and transform data contained in these files into
        structured nodes and relationships

        :param source: a given source or database listed in the staging area
        :param version: a given version of the selected database
        :param files: a dictionary/mapping-like object containing attribute_name-file_path pairs.
        Files parameter should map a given attribute to the respective file and pandas DataFrame.
        """
        if not source:
            source = ''

        if not version:
            version = ''

        self._source = source
        self._version = version
        self._df = pd.DataFrame()
        self._files = {}
        self._attrs = {}

        for key, file in files.items():

            if hasattr(self, key):
                raise ValueError(f'Invalid input for file {file} with key {key}')

            file_path = os.path.join(STAGING_AREA_PATH, source, version, file)

            if os.path.exists(file_path):

                df = pd.DataFrame()

                self._files[key] = file_path
                self.attrs = (key, df)

            else:
                raise FileNotFoundError(f'Could not found file {file_path} with key {key}')

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def df(self) -> pd.DataFrame:
        return self._df

    @property
    def files(self) -> Dict[str, str]:
        return self._files.copy()

    # --------------------------------------------------------
    # Dynamic properties
    # --------------------------------------------------------
    @property
    def attrs(self) -> Dict[str, pd.DataFrame]:
        return {key: df.copy() for key, df in self._attrs.items()}

    @attrs.setter
    def attrs(self, value: Tuple[str, pd.DataFrame]):

        key, df = value

        if self.is_valid_key(key) and isinstance(df, pd.DataFrame):
            self._attrs[key] = df
            setattr(self, key, df)
            return

        raise AttributeError(f'Attribute with name {key} cannot be set')

    # --------------------------------------------------------
    # Transformer Python API
    # --------------------------------------------------------
    def __str__(self):
        return self._df.__str__()

    def __getattr__(self, item) -> pd.DataFrame:
        if self.is_valid_key(item):
            df = pd.DataFrame()
            self.attrs = (item, df)

        return self._attrs[item]

    def __setattr__(self, key, value):
        if self.is_valid_key(key):
            self.attrs = (key, value)

    def __getitem__(self, item) -> pd.DataFrame:
        return self._attrs.__getitem__(item)

    def __setitem__(self, key, value):
        if self.is_valid_key(key):
            self.attrs = (key, value)

    def keys(self):
        return self._attrs.keys()

    def values(self):
        return self._attrs.values()

    def items(self):
        return self._attrs.items()

    def get(self, key: str, default=None) -> pd.DataFrame:
        return self._attrs.get(key, default)

    def update(self, values: Dict[str, pd.DataFrame]):

        for key, val in values.items():
            self.attrs = (key, val)

    # --------------------------------------------------------
    # Transformer API
    # --------------------------------------------------------
    @abstractmethod
    def read(self, *args, **kwargs):

        """
        The method responsible for reading files' mapping-like object into pandas DataFrames.
        The pandas DataFrames are loaded into the attrs attribute.
        Additionally, the transformer instance is loaded with the corresponding attributes

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer children.

        :param args:
        :param kwargs:
        :return:
        """

        pass

    @abstractmethod
    def process(self, *args, **kwargs):
        pass

    @abstractmethod
    def integrate(self, *args, **kwargs):
        pass

    # ----------------------------------------
    # Utilities methods
    # ----------------------------------------
    def read_json_lines(self):

        for key, file_path in self._files.items():
            df = read_json_lines(file_path)
            self.attrs = (key, df)

    def is_valid_key(self, key: str) -> bool:

        if key not in self._attrs and key not in self.__dict__:
            return True

        elif key in self._attrs and hasattr(self, key):
            return True

        return False

    def node_snapshot(self) -> pd.DataFrame:
        df = self.node.node_to_df()
        df.reset_index(inplace=True, drop=True)
        return df
