import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import List, Type, Dict, Tuple, Callable

import pandas as pd

from protrend.io import read_json_frame, write_json_frame, read
from protrend.model import BaseNode
from protrend.utils import Settings, DefaultProperty, SetList
from protrend.utils.processors import apply_processors


class AbstractConnector(metaclass=ABCMeta):
    """
    Connector interface.
    The following methods must be implemented to set up a connector for each relationship
    """

    # --------------------------------------------------------
    # Connector API
    # --------------------------------------------------------
    @abstractmethod
    def connect(self):
        """
        The method responsible for connecting an integrated pandas DataFrames to the remaining database.
        The connect method should set up multiple connections into several pandas DataFrames
        using either protrend identifiers only

        Interface implementation with unknown signature.
        Concrete implementations are available at the connector children.

        :return:
        """
        pass

    def write(self):
        """
        The method responsible for writing the connected nodes/relationships pandas DataFrame
        to be connected in the neo4j database using neomodel and Node class extension.

        Interface implementation with unknown signature.
        Concrete implementations are available at the connector or connector's sub-classes.

        """
        pass


CONNECTOR_STACK = dict(source='integrated_source.json',
                       effector='integrated_effector.json',
                       evidence='integrated_evidence.json',
                       gene='integrated_gene.json',
                       operon='integrated_operon.json',
                       organism='integrated_organism.json',
                       pathway='integrated_pathway.json',
                       publication='integrated_publication.json',
                       regulator='integrated_regulator.json',
                       rfam='integrated_regulatoryfamily.json',
                       rin='integrated_regulatoryinteraction.json',
                       tfbs='integrated_tfbs.json',
                       promoter_region='integrated_promoterregion.json',
                       motif='integrated_motif.json')


class Connector(AbstractConnector):
    """
    A connector is responsible for reading, processing, integrating and writing relationships files.

    A connector starts with data from the data lake and ends with structured relationships.
    """
    source = DefaultProperty()
    version = DefaultProperty()

    from_node = DefaultProperty()
    to_node = DefaultProperty()

    def __init_subclass__(cls, **kwargs):

        source = kwargs.get('source')
        cls.source.set_default(cls, source)

        version = kwargs.get('version')
        cls.version.set_default(cls, version)

        from_node = kwargs.pop('from_node', None)
        cls.from_node.set_default(cls, from_node)

        to_node = kwargs.pop('to_node', None)
        cls.to_node.set_default(cls, to_node)

        register = kwargs.pop('register', False)

        if register:
            from protrend.pipeline import Pipeline
            Pipeline.register_connector(cls, **kwargs)

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 from_node: Type[BaseNode] = None,
                 to_node: Type[BaseNode] = None):
        """
        The connector object uses results obtained during the transformation procedures
        for a given neomodel node entity.
        This object is responsible for reading several files and digest/process them into relationship instances.

        A given source/database contained in the data lake must be provided
        together with the files required for the connection to take place.

        A pandas DataFrame is the main engine to read, load, process and transform data contained in these files into
        structured relationships

        :type source: str
        :type version: str
        :type from_node: Type[BaseNode]
        :type to_node: Type[Node]

        :param source: The name of the data source in the data lake (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the data lake (e.g. 0.0.0, 0.0.1, etc)
        :param from_node: The source node type associated with this connector, and thus the source of the relation.
        Note that, it should be created only one connector for each node-node relationship
        :param to_node: The target node type associated with this connector, and thus the end of the relation.
        Note that, it should be created only one connector for each node-node relationship
        """
        self.source = source
        self.version = version
        self.from_node = from_node
        self.to_node = to_node

        self._write_stack = []

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def write_path(self) -> str:
        return os.path.join(Settings.data_lake, self.source, self.version)

    # --------------------------------------------------------
    # Python API
    # --------------------------------------------------------
    def __str__(self):
        return f'{self.__class__.__name__}: {self.source} - {self.version} - ' \
               f'{self.from_node.node_name()} - {self.to_node.node_name()}'

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other: 'Connector'):

        other_values = {other.__class__.__name__, other.source, other.version, other.from_node.node_name(),
                        other.to_node.node_name()}

        if self.__class__.__name__ not in other_values:
            return False

        if self.source not in other_values:
            return False

        if self.version not in other_values:
            return False

        if self.from_node.node_name() not in other_values:
            return False

        if self.to_node.node_name() not in other_values:
            return False

        return True

    # --------------------------------------------------------
    # Connector API
    # --------------------------------------------------------
    @abstractmethod
    def connect(self):
        """
        The method responsible for connecting an integrated pandas DataFrames to the remaining database.
        The connect method should set up multiple connections into several pandas DataFrames
        using either protrend identifiers only

        Interface implementation with unknown signature.
        Concrete implementations are available at the connector children.

        :return:
        """
        pass

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
    def connected_file(cls):
        from_node_name = cls.from_node.get_default(cls).node_name()
        to_node_name = cls.to_node.get_default(cls).node_name()
        return f'connected_{from_node_name}_{to_node_name}.json'

    def stack_connections(self, df: pd.DataFrame):
        df = df.copy()
        df = df.reset_index(drop=True)
        file_name = f'connected_{self.from_node.node_name()}_{self.to_node.node_name()}.json'
        fp = os.path.join(self.write_path, file_name)
        json_partial = partial(write_json_frame, file_path=fp, df=df)
        self._write_stack.append(json_partial)

    def transform_stacks(self,
                         source: str,
                         target: str,
                         source_column: str,
                         target_column: str,
                         source_on: str = None,
                         target_on: str = None,
                         source_processors: Dict[str, List[Callable]] = None,
                         target_processors: Dict[str, List[Callable]] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:

        default_source_cols = SetList([source_column])
        default_target_cols = SetList([target_column])

        if source_on:
            default_source_cols.append(source_on)

        if target_on:
            default_target_cols.append(target_on)

        if not source_processors:
            source_processors = {}

        if not target_processors:
            target_processors = {}

        default_source_cols.extend(list(source_processors.keys()))
        default_target_cols.extend(list(target_processors.keys()))

        file = CONNECTOR_STACK[source]
        source_df = read(source=self.source, version=self.version, file=file,
                         reader=read_json_frame, default=pd.DataFrame(columns=default_source_cols))
        source_df = apply_processors(source_df, **source_processors)
        source_df = source_df.rename(columns={source_column: 'source_col'})

        file = CONNECTOR_STACK[target]
        target_df = read(source=self.source, version=self.version, file=file,
                         reader=read_json_frame, default=pd.DataFrame(columns=default_target_cols))
        target_df = apply_processors(target_df, **target_processors)
        target_df = target_df.rename(columns={target_column: 'target_col'})

        return source_df, target_df

    @staticmethod
    def merge_source_target(source_df: pd.DataFrame,
                            target_df: pd.DataFrame,
                            source_on: str = None,
                            target_on: str = None,
                            cardinality: str = 'many_to_many') -> Tuple[List[str], List[str]]:
        if source_on:
            df = pd.merge(source_df, target_df, left_on=source_on, right_on=target_on)
            df = df.drop_duplicates(subset=['source_col', 'target_col'])

        else:
            source_df = source_df.reset_index(drop=True)
            target_df = target_df.reset_index(drop=True)
            df = pd.concat([source_df, target_df], axis=1)
            df = df.drop_duplicates(subset=['source_col', 'target_col'])

        if cardinality == 'many_to_many':
            source_ids = df['source_col'].to_list()
            target_ids = df['target_col'].to_list()

        elif cardinality == 'one_to_many':
            target_ids = df['target_col'].to_list()

            source_ids = df['source_col'].dropna().to_list()
            source_ids *= len(target_ids)

        elif cardinality == 'many_to_one':
            source_ids = df['source_col'].to_list()

            target_ids = df['target_col'].dropna().to_list()
            target_ids *= len(source_ids)

        else:
            source_ids = []
            target_ids = []

        return source_ids, target_ids

    def connection_frame(self,
                         source_ids: List[str],
                         target_ids: List[str],
                         kwargs: Dict[str, List[str]] = None) -> pd.DataFrame:

        size = len(source_ids)

        connection = {'from_node': [self.from_node.node_name()] * size,
                      'to_node': [self.to_node.node_name()] * size,
                      'from_identifier': source_ids.copy(),
                      'to_identifier': target_ids.copy(),
                      'load': ['create'] * size,
                      'what': ['relationships'] * size}

        if not kwargs:
            kwargs = {}

        connection.update(kwargs)

        return pd.DataFrame(connection)

    def create_connection(self,
                          source: str,
                          target: str,
                          source_column: str = 'protrend_id',
                          target_column: str = 'protrend_id',
                          cardinality: str = 'many_to_many',
                          source_on: str = None,
                          target_on: str = None,
                          source_processors: Dict[str, List[Callable]] = None,
                          target_processors: Dict[str, List[Callable]] = None) -> pd.DataFrame:

        source_df, target_df = self.transform_stacks(source=source,
                                                     target=target,
                                                     source_column=source_column,
                                                     target_column=target_column,
                                                     source_on=source_on,
                                                     target_on=target_on,
                                                     source_processors=source_processors,
                                                     target_processors=target_processors)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          source_on=source_on, target_on=target_on,
                                                          cardinality=cardinality)

        return self.connection_frame(source_ids=source_ids, target_ids=target_ids)
