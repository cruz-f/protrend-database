import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import Dict, Tuple, Union

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
        self._files = {}
        self._attrs = {}
        self._write_stack = []

        for key, file in files.items():

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

        self._attrs[key] = df
        return

    @property
    def write_path(self):
        return os.path.join(DATA_LAKE_PATH, self._source, self._version)

    # --------------------------------------------------------
    # Transformer Python API
    # --------------------------------------------------------
    def __str__(self):
        return self.attrs.__str__()

    def __getitem__(self, item) -> pd.DataFrame:
        return self._attrs.__getitem__(item)

    def __setitem__(self, key, value):
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
    def transform(self, *args, **kwargs):
        pass

    @abstractmethod
    def load_nodes(self, *args, **kwargs):
        pass

    @abstractmethod
    def load_relationships(self, *args, **kwargs):
        # relationship dataframe structure
        # from
        # to
        # from_property
        # to_property

        # relationship dataframe custom structure
        # name
        # url
        # external_identifier
        pass

    def write(self):

        if not os.path.exists(self.write_path):
            os.makedirs(self.write_path)

        for csv in self._write_stack:
            csv()

        self._write_stack = []

    # ----------------------------------------
    # Utilities methods
    # ----------------------------------------
    def read_json_lines(self):

        for key, file_path in self._files.items():
            df = read_json_lines(file_path)
            self.attrs = (key, df)

    def stack_csv(self, name: str, df: pd.DataFrame):

        df_copy = df.copy(deep=True)
        fp = os.path.join(self.write_path, f'{name}.csv')
        csv = partial(df_copy.to_csv, path_or_buf=fp)
        self._write_stack.append(csv)

    def last_node(self) -> Union['Node', None]:
        return self.node.last_node()

    def node_snapshot(self) -> pd.DataFrame:
        df = self.node.node_to_df()
        df.reset_index(inplace=True, drop=True)
        return df
