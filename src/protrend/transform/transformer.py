import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import List, Type, Callable, Dict, Tuple

import neo4j.exceptions
import pandas as pd

from protrend.io import write_json_frame
from protrend.model import BaseNode, protrend_id_decoder, protrend_id_encoder
from protrend.utils import Settings, DefaultProperty, SetList
from protrend.utils.processors import apply_processors
from .transformations import drop_empty_string, drop_duplicates


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
    source = DefaultProperty()
    version = DefaultProperty()
    node = DefaultProperty()
    order = DefaultProperty()

    columns = SetList()

    def __init_subclass__(cls, **kwargs):
        source = kwargs.get('source')
        cls.source.set_default(cls, source)

        version = kwargs.get('version')
        cls.version.set_default(cls, version)

        node = kwargs.get('node')
        cls.node.set_default(cls, node)

        order = kwargs.get('order')
        cls.order.set_default(cls, order)

        register = kwargs.get('register', True)

        if register:
            from protrend.pipeline import Pipeline
            Pipeline.register_transformer(transformer=cls, source=source, version=version)

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 node: Type[BaseNode] = None,
                 order: int = None):

        """
        The transform object must implement several transformation procedures for a given neomodel Node.
        Thus, the transform must be able to read staging area files and digest/process/transform/integrate/write
        them into fully structured Node instances.

        A TransformerSettings object must be provided to a transformer.
        These settings set up a given data source contained in the staging area for the transformation to take place.

        Pandas is the main engine to read, process, transform and integrate data
        contained in the staging area files into structured nodes.

        :type source: str
        :type version: str
        :type node: Type[BaseNode]
        :type order: int

        :param source: The name of the data source in the staging area (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the staging area (e.g. 0.0.0, 0.0.1, etc)
        :param node: The node type associated with this transformer.
        Note that, it should be created only one transformer for each node
        :param order: The order value is used by the director to rearrange transformers and their execution.
        Transformers having higher order are executed first.
        """

        self.source = source
        self.version = version
        self.node = node
        self.order = order

        self._write_stack = []

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def node_factors(self) -> Dict[str, List[Callable]]:
        return self.node.node_factors

    @property
    def write_path(self) -> str:
        return os.path.join(Settings.data_lake, self.source, self.version)

    # --------------------------------------------------------
    # Python API
    # --------------------------------------------------------
    def __str__(self):
        return f'{self.__class__.__name__}: {self.source} - {self.version} - ' \
               f'{self.node.node_name()} - {list(self.node_factors.keys())}'

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

    def integrate(self, df: pd.DataFrame):

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
        # take a db snapshot for the current node and ensure uniqueness
        view = self.node_view()
        view, _ = self.factors_normalization(df=view)

        # ensure uniqueness
        df, factorized_cols = self.factors_normalization(df=df)

        # assign the integration and load columns
        df = df.assign(protrend_id=None, load=None, what='nodes')

        # try to integrate the new nodes by the node factors
        for factor in factorized_cols:
            mapper = view.dropna(subset=[factor])
            mapper = mapper.set_index(factor)
            mapper = mapper['protrend_id']
            to_fill = df[factor].map(mapper)
            filled = df['protrend_id'].fillna(to_fill)
            df = df.assign(protrend_id=filled)

        # those that mismatch the integration by the node factors (that is, the protrend id is missing) will be assigned
        # for creation
        mask = df['protrend_id'].isna()
        df.loc[mask, 'load'] = 'create'
        df.loc[~mask, 'load'] = 'update'

        # inferring the batch size, creating the new ids and assign to the df
        nodes = view['protrend_id'].to_list()
        last_node_idx = self.last_node_index(nodes)
        batch_size = mask.sum()

        batch_ids = self.protrend_identifiers_batch(last_node_idx=last_node_idx, size=batch_size)
        df.loc[mask, 'protrend_id'] = batch_ids

        df = df.drop(columns=factorized_cols)

        self.stack_integrated_nodes(df)
        self.stack_nodes(df)

    def write(self):

        if not os.path.exists(self.write_path):
            os.makedirs(self.write_path)

        for json_partial in self._write_stack:
            json_partial()

        self._write_stack = []

    # ----------------------------------------
    # Utilities
    # ----------------------------------------
    @classmethod
    def nodes_file(cls):
        node_name = cls.node.get_default(cls).node_name()
        return f'nodes_{node_name}.json'

    def empty_frame(self) -> pd.DataFrame:
        cols = [col for col in self.columns if col != 'protrend_id']
        return pd.DataFrame(columns=cols)

    def integrated_empty_frame(self) -> pd.DataFrame:
        cols = [col for col in self.columns]
        cols += ['load', 'what', 'node']
        return pd.DataFrame(columns=cols)

    def stack_json(self, name: str, df: pd.DataFrame):
        df = df.copy()
        df = df.reset_index(drop=True)
        fp = os.path.join(self.write_path, f'{name}.json')
        json_partial = partial(write_json_frame, file_path=fp, df=df)
        self._write_stack.append(json_partial)

    def stack_nodes(self, df: pd.DataFrame):
        df = df.copy()
        df = df.reset_index(drop=True)

        if df.empty:
            df = self.integrated_empty_frame()

        else:
            df.loc[:, 'node'] = self.node.node_name()

        node_cols = list(self.node.node_keys()) + ['load', 'what', 'node']
        cols_to_drop = [col for col in df.columns if col not in node_cols]
        df = df.drop(columns=cols_to_drop)

        df_name = f'nodes_{self.node.node_name()}'
        self.stack_json(df_name, df)

    def stack_integrated_nodes(self, df: pd.DataFrame):
        df = df.copy()
        df = df.reset_index(drop=True)

        if df.empty:
            df = self.integrated_empty_frame()

        else:
            df.loc[:, 'node'] = self.node.node_name()

        df_name = f'integrated_{self.node.node_name()}'
        self.stack_json(df_name, df)

    def stack_transformed_nodes(self, df: pd.DataFrame):
        df = df.copy()
        df = df.reset_index(drop=True)

        if df.empty:
            df = self.empty_frame()
        df_name = f'transformed_{self.node.node_name()}'
        self.stack_json(df_name, df)

    def factors_normalization(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, List[str]]:
        # renaming the factors
        to_assign = {}
        to_process = {}
        new_cols = []

        for i, (factor, processing_fns) in enumerate(self.node_factors.items()):
            new_col = f'node_factor_{i}'

            to_assign[new_col] = df[factor].copy()
            to_process[new_col] = processing_fns.copy()
            new_cols.append(new_col)

        df = df.assign(**to_assign)
        df = apply_processors(df, **to_process)
        df = drop_empty_string(df, *new_cols)
        df = drop_duplicates(df=df, subset=new_cols)
        return df, new_cols

    def protrend_identifiers_batch(self, last_node_idx: int, size: int) -> List[str]:

        return [protrend_id_encoder(self.node.header, self.node.entity, i)
                for i in range(last_node_idx + 1, last_node_idx + size + 1)]

    @staticmethod
    def last_node_index(nodes: List[str]) -> int:
        if not nodes:
            return 0

        sorted_nodes = sorted(nodes, key=protrend_id_decoder, reverse=True)
        return protrend_id_decoder(sorted_nodes[0])

    def node_view(self) -> pd.DataFrame:
        try:
            return self.node.node_to_df()

        except (neo4j.exceptions.Neo4jError, neo4j.exceptions.DriverError):
            return pd.DataFrame(columns=list(self.node.node_keys()))
