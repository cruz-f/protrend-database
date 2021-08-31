import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import Tuple, Union, List, Type, Callable, Dict, Sequence, Set

import pandas as pd

from protrend.model.node import Node, protrend_id_decoder, protrend_id_encoder
from protrend.transform.processors import take_last
from protrend.transform.settings import TransformerSettings
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

    default_settings: Type[TransformerSettings] = TransformerSettings
    columns: Set[str] = set()

    def __init__(self, settings: TransformerSettings = None):

        """
        The transform object must implement several transformation procedures for a given neomodel Node.
        Thus, the transform must be able to read staging area files and digest/process/transform/integrate/write
        them into fully structured Node instances.

        A TransformerSettings object must be provided to a transformer.
        These settings set up a given data source contained in the staging area for the transformation to take place.

        Pandas is the main engine to read, process, transform and integrate data
        contained in the staging area files into structured nodes.

        :type settings: TransformerSettings
        :param settings: a TransformerSettings than contains all settings for source,
        version and files to perform the transformation on
        """
        if not settings:
            settings = self.default_settings()

        self._settings = settings

        self._transform_stack = {}
        self._write_stack = []

        self._load_settings()

    def _load_settings(self):

        for key, file in self._settings.transform.items():

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
    def node(self) -> Type[Node]:
        return self.settings.node

    @property
    def node_factors(self) -> Tuple[str]:
        return self.settings.node_factors

    @property
    def transform_stack(self) -> Dict[str, str]:
        return self._transform_stack

    @property
    def write_stack(self) -> List[Callable]:
        return self._write_stack

    @property
    def order(self) -> int:
        return self.settings.order

    @property
    def write_path(self) -> str:
        return os.path.join(DATA_LAKE_PATH, self.source, self.version)

    # --------------------------------------------------------
    # Python API
    # --------------------------------------------------------
    def __str__(self):
        return f'{self.__class__.__name__}: {self.settings}'

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other: 'Transformer'):

        return self.__class__.__name__ and other.__class__.__name__ and self.settings == other.settings

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

        # find/set protrend identifiers for update nodes
        ids_mask = self.find_snapshot(nodes=nodes, snapshot=snapshot, node_factors=self.node_factors)
        nodes['protrend_id'] = snapshot.loc[ids_mask, 'protrend_id']
        nodes['load'] = 'update'
        nodes['what'] = 'nodes'

        return nodes

    def _create_nodes(self, df: pd.DataFrame, mask: pd.Series) -> pd.DataFrame:

        # nodes to be created
        nodes = df[~mask]

        # create/set new protrend identifiers
        create_size, _ = nodes.shape
        create_identifiers = self.protrend_identifiers_batch(create_size)
        nodes['protrend_id'] = create_identifiers
        nodes['load'] = 'create'
        nodes['what'] = 'nodes'

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
        df = self.drop_duplicates(df=df, subset=self.node_factors, perfect_match=False, preserve_nan=True)
        df = df.reset_index()

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

        for csv in self._write_stack:
            csv()

        self._write_stack = []

    # ----------------------------------------
    # Utilities methods
    # ----------------------------------------
    def empty_frame(self) -> pd.DataFrame:
        cols = [col for col in self.columns if col != 'protrend_id']
        return pd.DataFrame(columns=cols)

    def stack_csv(self, name: str, df: pd.DataFrame):
        df_copy = df.copy(deep=True)
        fp = os.path.join(self.write_path, f'{name}.csv')
        csv = partial(df_copy.to_csv, path_or_buf=fp, index=False)
        self._write_stack.append(csv)

    def _stack_nodes(self, df: pd.DataFrame):
        node_cols = list(self.node.node_keys())
        cols_to_drop = [col for col in df.columns if col not in node_cols]
        df = df.drop(columns=cols_to_drop)
        df_name = f'nodes_{self.node.node_name()}'
        self.stack_csv(df_name, df)

    def _stack_integrated_nodes(self, df: pd.DataFrame):
        df_name = f'integrated_{self.node.node_name()}'
        self.stack_csv(df_name, df)

    def _stack_transformed_nodes(self, df: pd.DataFrame):
        if df.empty:
            df = self.empty_frame()
        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

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

        return df

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
                for i in range(integer + 1, size + 1)]

    def last_node(self) -> Union['Node', None]:
        return self.node.last_node()

    def node_view(self) -> pd.DataFrame:
        df = self.node.node_to_df()
        return df
