import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import Tuple, Union, List, Type, Callable, Dict, Set

import pandas as pd

from protrend.io.json import write_json_frame
from protrend.model.node import Node, protrend_id_decoder, protrend_id_encoder
from protrend.transform.processors import take_last, apply_processors, to_nan
from protrend.utils.miscellaneous import is_null
from protrend.utils.settings import STAGING_AREA_PATH, DATA_LAKE_PATH


class AbstractTransformer(metaclass=ABCMeta):
    """
    Transformer interface.
    The following methods must be implemented to set up a transformer for each node
    """

    # --------------------------------------------------------
    # Transformer API
    # --------------------------------------------------------
    @abstractmethod
    def transform(self) -> pd.DataFrame:

        """
        The method responsible for transforming multiple pandas DataFrames into an annotated, cleaned and standardized
        pandas DataFrame ready to be integrated into structured nodes.

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer or transformer's sub-classes.

        :return:
        """

        pass

    @abstractmethod
    def integrate(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        The method responsible for integrating the transformed pandas DataFrame with the current state of the neo4j
        database using neomodel and Node class extension.

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer or transformer's sub-classes.

        :param df: transformed pandas DataFrame
        :return: it creates a new pandas DataFrame of the integrated data
        """
        pass

    def write(self):
        """
        The method responsible for writing the transformed, integrated and nodes pandas DataFrame
        to be loaded and connected in the neo4j database using neomodel and Node class extension.

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer or transformer's sub-classes.

        """
        pass


class Transformer(AbstractTransformer):
    """
    A transformer is responsible for reading, processing, transforming, integrating and writing node files.
    A transformer starts with data from the staging area and ends with structured nodes.

    This Transformer object has implemented several utilities for the hard-working de facto transformers to be created
    for each data source and node
    """
    default_transform_stack: Dict[str, str] = {}
    default_source: str = ''
    default_version: str = '0.0.0'
    default_node: Type[Node] = Node
    default_node_factors: Tuple[str] = ()
    default_order: int = 0
    columns: Set[str] = set()
    read_columns: Set[str] = set()

    def __init__(self,
                 transform_stack: Dict[str, str] = None,
                 source: str = None,
                 version: str = None,
                 node: Type[Node] = None,
                 node_factors: Tuple[str] = None,
                 order: int = None):

        """
        The transform object must implement several transformation procedures for a given neomodel Node.
        Thus, the transform must be able to read staging area files and digest/process/transform/integrate/write
        them into fully structured Node instances.

        A TransformerSettings object must be provided to a transformer.
        These settings set up a given data source contained in the staging area for the transformation to take place.

        Pandas is the main engine to read, process, transform and integrate data
        contained in the staging area files into structured nodes.

        :type transform_stack: Dict[str, str]
        :type source: str
        :type version: str
        :type node: Type[Node]
        :type node_factors: Tuple[str]
        :type order: int

        :param transform_stack: Dictionary containing the pair name and file name.
        The key should be used to identify the file in the transform stack,
        whereas the value should be the file name in the staging area or data lake
        :param source: The name of the data source in the staging area (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the staging area (e.g. 0.0.0, 0.0.1, etc)
        :param node: The node type associated with this transformer.
        Note that, it should be created only one transformer for each node
        :param node_factors: The node attributes that must be used during the integration.
        The node factors will be used iteratively one-by-one to find nodes already available in the database.
        For instance, if the Regulator uniprot_accession attribute is used as node factor,
        regulators in transformation and available in the database will be merged by their uniprot_accession.
        :param order: The order value is used by the director to rearrange transformers and their execution.
        Transformers having higher order are executed first.
        """

        self._transform_stack = {}
        self._write_stack = []
        self._source = source
        self._version = version
        self._node = node
        self._node_factors = node_factors
        self._order = order

        self.load_transform_stack(transform_stack)

    def load_transform_stack(self, transform_stack: Dict[str, str] = None):

        self._transform_stack = {}

        if not transform_stack:
            transform_stack = self.default_transform_stack

        for key, file in transform_stack.items():

            sa_file = os.path.join(STAGING_AREA_PATH, self.source, self.version, file)
            dl_file = os.path.join(DATA_LAKE_PATH, self.source, self.version, file)

            if os.path.exists(sa_file):

                self._transform_stack[key] = sa_file

            else:

                self._transform_stack[key] = dl_file

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def source(self) -> str:
        if not self._source:
            return self.default_source

        return self._source

    @property
    def version(self) -> str:
        if not self._version:
            return self.default_version

        return self._version

    @property
    def node(self) -> Type[Node]:
        if not self._node:
            return self.default_node

        return self._node

    @property
    def node_factors(self) -> Tuple[str]:
        if not self._node_factors:
            return self.default_node_factors

        return self._node_factors

    @property
    def transform_stack(self) -> Dict[str, str]:
        return self._transform_stack

    @property
    def order(self) -> int:
        if self._order is None:
            return self.default_order

        return self._order

    @property
    def write_path(self) -> str:
        return os.path.join(DATA_LAKE_PATH, self.source, self.version)

    # --------------------------------------------------------
    # Python API
    # --------------------------------------------------------
    def __str__(self):
        return f'{self.__class__.__name__}: {self.source} - {self.version} - ' \
               f'{self.node.node_name()} - {self.node_factors}'

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other: 'Transformer'):

        other_values = {other.__class__.__name__, other.source, other.version, other.node.node_name(),
                        other.node_factors}

        if self.__class__.__name__ not in other_values:
            return False

        if self.source not in other_values:
            return False

        if self.version not in other_values:
            return False

        if self.node.node_name() not in other_values:
            return False

        if self.node_factors not in other_values:
            return False

        return True

    # --------------------------------------------------------
    # Transformer API
    # --------------------------------------------------------
    @abstractmethod
    def transform(self):

        """
        The method responsible for transforming multiple pandas DataFrames into an annotated, cleaned and standardized
        pandas DataFrame ready to be integrated and transformed into structured nodes.

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer children.

        :return:
        """

        pass

    def _update_nodes(self, df: pd.DataFrame, mask: pd.Series, snapshot: pd.DataFrame) -> pd.DataFrame:

        # nodes to be updated
        nodes = df[mask]

        if nodes.empty:
            nodes['protrend_id'] = None
            nodes['load'] = None
            nodes['what'] = None
            return nodes

        # find/set protrend identifiers for update nodes
        ids_mask = self.find_snapshot(nodes=nodes, snapshot=snapshot, node_factors=self.node_factors)
        nodes.loc[:, 'protrend_id'] = snapshot.loc[ids_mask, 'protrend_id']
        nodes.loc[:, 'load'] = 'update'
        nodes.loc[:, 'what'] = 'nodes'

        return nodes

    def _create_nodes(self, df: pd.DataFrame, mask: pd.Series) -> pd.DataFrame:

        # nodes to be created
        nodes = df[~mask]
        if nodes.empty:
            nodes['protrend_id'] = None
            nodes['load'] = None
            nodes['what'] = None
            return nodes

        # create/set new protrend identifiers
        create_size, _ = nodes.shape
        create_identifiers = self.protrend_identifiers_batch(create_size)
        nodes.loc[:, 'protrend_id'] = create_identifiers
        nodes.loc[:, 'load'] = 'create'
        nodes.loc[:, 'what'] = 'nodes'

        return nodes

    def integrate(self, df: pd.DataFrame) -> pd.DataFrame:

        """
        The method responsible for integrating the transformed pandas DataFrame with the current state of the neo4j
        database using neomodel and node class extension

        The node factors are node attributes used for data normalization and integration.
        The node factors are used to retrieve the protrend identifier associated with each row of the database.
        Pandas boolean masks are used to query the database improving search performance.

        The nodes identified by the distinct node factors are marked as update nodes. Thus, these nodes identified by
        their protrend identifier will be updated with new data.

        The nodes missing a hit with the current state of the database will be created according
        to the data contained in the DataFrame.

        :param df: transformed pandas DataFrame
        :return: it creates a new pandas DataFrame of the integrated data
        """
        # ensure uniqueness
        df = self.standardize_nulls(df=df)
        df = self.drop_duplicates(df=df, subset=self.node_factors, perfect_match=False, preserve_nan=True)

        # take a db snapshot for the current node
        snapshot = self.node_view()

        # find matching nodes according to several node factors/properties
        mask = self.find_nodes(nodes=df, snapshot=snapshot, node_factors=self.node_factors)

        # nodes to be updated
        update_nodes = self._update_nodes(df=df, mask=mask, snapshot=snapshot)

        # nodes to be created
        create_nodes = self._create_nodes(df=df, mask=mask)

        # concat both dataframes
        df = pd.concat([create_nodes, update_nodes])

        self._stack_integrated_nodes(df)
        self._stack_nodes(df)

        return df

    def write(self):

        if not os.path.exists(self.write_path):
            os.makedirs(self.write_path)

        for json_partial in self._write_stack:
            json_partial()

        self._write_stack = []

    # ----------------------------------------
    # Utilities
    # ----------------------------------------
    def empty_frame(self) -> pd.DataFrame:
        cols = [col for col in self.columns if col != 'protrend_id']
        return pd.DataFrame(columns=cols)

    def stack_json(self, name: str, df: pd.DataFrame):
        df = df.copy(deep=True)
        df = df.reset_index(drop=True)
        fp = os.path.join(self.write_path, f'{name}.json')
        json_partial = partial(write_json_frame, file_path=fp, df=df)
        self._write_stack.append(json_partial)

    def _stack_nodes(self, df: pd.DataFrame):
        if df.empty:
            df = self.empty_frame()
            df.loc[:, 'load'] = None
            df.loc[:, 'what'] = None

        df.loc[:, 'node'] = self.node.node_name()

        node_cols = list(self.node.node_keys()) + ['load', 'what', 'node']
        cols_to_drop = [col for col in df.columns if col not in node_cols]
        df = df.drop(columns=cols_to_drop)

        df_name = f'nodes_{self.node.node_name()}'
        self.stack_json(df_name, df)

    def _stack_integrated_nodes(self, df: pd.DataFrame):
        if df.empty:
            df = self.empty_frame()
            df.loc[:, 'load'] = None
            df.loc[:, 'what'] = None

        df.loc[:, 'node'] = self.node.node_name()

        df_name = f'integrated_{self.node.node_name()}'
        self.stack_json(df_name, df)

    def _stack_transformed_nodes(self, df: pd.DataFrame):
        if df.empty:
            df = self.empty_frame()
        df_name = f'transformed_{self.node.node_name()}'
        self.stack_json(df_name, df)

    def standardize_nulls(self, df: pd.DataFrame) -> pd.DataFrame:

        processors = {col: to_nan for col in self.node_factors}

        df = apply_processors(df, **processors)
        return df

    @staticmethod
    def create_input_value(df: pd.DataFrame, col: str) -> pd.DataFrame:
        df['input_value'] = df[col]
        return df

    @staticmethod
    def select_columns(df: pd.DataFrame, *columns: str) -> pd.DataFrame:
        df = df[list(columns)]
        return df

    @staticmethod
    def merge_columns(df: pd.DataFrame, column: str, left: str, right: str) -> pd.DataFrame:
        df[column] = df[left].fillna(df[right])
        df = df.drop(columns=[left, right])
        return df

    @staticmethod
    def concat_columns(df: pd.DataFrame, column: str, left: str, right: str) -> pd.DataFrame:
        df[column] = df[left] + df[right]
        df = df.drop(columns=[left, right])
        return df

    @staticmethod
    def group_by(df: pd.DataFrame,
                 column: str,
                 aggregation: Dict[str, Callable],
                 default: Callable = take_last) -> pd.DataFrame:

        agg = {col: default for col in df.columns if col != column}
        agg.update(aggregation)

        df = df.groupby(df[column]).aggregate(agg)
        df = df.reset_index()
        return df

    @staticmethod
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

    @staticmethod
    def find_nodes(nodes: pd.DataFrame, snapshot: pd.DataFrame, node_factors: Tuple[str]) -> pd.Series:

        n_rows, _ = nodes.shape

        mask = pd.Series([False] * n_rows)

        factors_masks = []
        for factor in node_factors:

            if factor in nodes.columns and factor in snapshot.columns:
                snapshot_values = snapshot[factor]
                is_in_mask = nodes[factor].isin(snapshot_values)
                is_not_null_mask = ~ nodes[factor].map(is_null)
                factor_mask = (is_in_mask & is_not_null_mask)
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
                is_in_mask = snapshot[factor].isin(node_values)
                is_not_null_mask = ~ snapshot[factor].map(is_null)
                factor_mask = (is_in_mask & is_not_null_mask)
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
                for i in range(integer + 1, integer + size + 1)]

    def last_node(self) -> Union['Node', None]:
        return self.node.last_node()

    def node_view(self) -> pd.DataFrame:
        df = self.node.node_to_df()
        return df
