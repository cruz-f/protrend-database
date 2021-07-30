import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import Dict, Tuple, Union, Callable, List

import pandas as pd

from protrend.io.json import read_json_lines
from protrend.model.node import Node, protrend_id_decoder, protrend_id_encoder
from protrend.transform.settings import TransformerSettings
from protrend.utils.settings import STAGING_AREA_PATH, DATA_LAKE_PATH


class Transformer(metaclass=ABCMeta):
    """
    Transformer interface.

    The following methods must be implemented to set up a transformer for each node
    and relationship

    A transformer is responsible for reading, processing, integrating and writing node and relationships files.

    A transformer starts with data from the staging area and ends with structured nodes and relationships.
    """

    def __init__(self, settings: TransformerSettings):

        """
        The transform object contains and uses transformation procedures for a given neomodel node entity.
        This object is responsible for reading several files and digest/process them into node instances.

        A given source/database contained in the staging area must be provided
        together with the files required for the transformation to take place.

        An attribute will be created for each key-value pair in the files' mapping-like object.

        A pandas DataFrame is the main engine to read, load, process and transform data contained in these files into
        structured nodes and relationships

        :param settings: a TransformerSettings than contains all settings for source,
        version and files to perform the transformation on
        """

        self._settings = settings

        self._files = {}
        self._attrs = {}
        self._write_stack = []

        self._load_settings()

    def _load_settings(self):

        for key, file in self._settings.files.items():

            file_path = os.path.join(STAGING_AREA_PATH, self.source, self.version, file)

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
    def settings(self) -> TransformerSettings:
        return self._settings

    # --------------------------------------------------------
    # Dynamic properties
    # --------------------------------------------------------
    @property
    def source(self) -> str:
        return self.settings.source

    @property
    def version(self) -> str:
        return self.settings.version

    @property
    def files(self) -> Dict[str, str]:
        return self._files

    @property
    def attrs(self) -> Dict[str, pd.DataFrame]:
        return self._attrs

    @attrs.setter
    def attrs(self, value: Tuple[str, pd.DataFrame]):

        key, df = value

        self._attrs[key] = df
        return

    @property
    def node(self) -> Node:
        return self.settings.node

    @property
    def node_factors(self) -> Node:
        return self.settings.node_factors

    @property
    def write_path(self) -> str:
        return os.path.join(DATA_LAKE_PATH, self.source, self.version)

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

    def _load_nodes(self, df: pd.DataFrame, identifiers: List[str], node_factory: Callable) -> pd.DataFrame:

        df[self.node.identifying_property] = identifiers
        node_factory(nodes=df, save=True)

        return df

    def load_nodes(self, df: pd.DataFrame, *node_factors: Tuple[str]) -> pd.DataFrame:

        if not node_factors:
            node_factors = self.node_factors

        # take a db snapshot for the current node
        snapshot = self.node_snapshot()

        # find matching nodes according to several node factors/properties
        nodes_mask = self.find_nodes(nodes=df, snapshot=snapshot, node_factors=node_factors)

        # nodes to be updated
        update_nodes = df[nodes_mask]

        snapshot_mask = self.find_snapshot(nodes=update_nodes, snapshot=snapshot, node_factors=node_factors)
        update_identifiers = snapshot.loc[snapshot_mask, self.node.identifying_property]

        update_nodes = self._load_nodes(df=update_nodes,
                                        identifiers=update_identifiers,
                                        node_factory=self.node.node_update_from_df)

        # nodes to be created
        create_nodes = df[~nodes_mask]

        size, _ = create_nodes.shape
        create_identifiers = self.protrend_identifiers_batch(size)

        create_nodes = self._load_nodes(df=create_nodes, identifiers=create_identifiers,
                                        node_factory=self.node.node_from_df)

        df = pd.concat([create_nodes, update_nodes], axis=0)
        self.stack_csv(self.node.node_name(), df)

        return df

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

    @staticmethod
    def find_nodes(nodes: pd.DataFrame, snapshot: pd.DataFrame, node_factors: Tuple[str]) -> pd.Series:

        n_rows, _ = nodes.shape

        mask = pd.Series([False] * n_rows)

        factors_masks = []
        for factor in node_factors:

            if factor in nodes.columns and factor in snapshot.columns:
                snapshot_values = snapshot[factor]
                factor_mask = nodes[factor].isin(snapshot_values)
                factors_masks.append(factor_mask)

        if factors_masks:
            mask = pd.concat(factors_masks, axis=1)
            mask = mask.any(axis=1)

        return mask

    @staticmethod
    def find_snapshot(nodes: pd.DataFrame, snapshot: pd.DataFrame, node_factors: Tuple[str]) -> pd.Series:

        n_rows, _ = snapshot.shape

        mask = pd.Series([False] * n_rows)

        factors_masks = []
        for factor in node_factors:

            if factor in nodes.columns and factor in snapshot.columns:
                node_values = nodes[factor]
                factor_mask = snapshot[factor].isin(node_values)
                factors_masks.append(factor_mask)

        if factors_masks:
            mask = pd.concat(factors_masks, axis=1)
            mask = mask.any(axis=1)

        return mask

    def protrend_identifiers_batch(self, size):

        last_node = self.node.last_node()

        if last_node is None:
            integer = 0
        else:
            integer = protrend_id_decoder(last_node.protrend_id)

        return [protrend_id_encoder(self.node.header, self.node.entity, i)
                for i in range(integer+1, size+1)]

    @staticmethod
    def make_relationship_df(n_rows: int,
                             from_node: str,
                             to_node: str,
                             from_property: str,
                             to_property: str) -> pd.DataFrame:

        relationship = {'from': [from_node] * n_rows,
                        'to': [to_node] * n_rows,
                        'from_property': [from_property] * n_rows,
                        'to_property': [to_property] * n_rows}

        return pd.DataFrame(relationship)

    def last_node(self) -> Union['Node', None]:
        return self.node.last_node()

    def node_snapshot(self) -> pd.DataFrame:
        df = self.node.node_to_df()
        return df
