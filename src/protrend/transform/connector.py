import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import List, Any, Type, Dict

import pandas as pd

from protrend.io.json import write_json_frame
from protrend.model.node import Node
from protrend.utils import Settings, WriteStack, DefaultProperty


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


class Connector(AbstractConnector):
    """
    A connector is responsible for reading, processing, integrating and writing relationships files.

    A connector starts with data from the data lake and ends with structured relationships.
    """
    source = DefaultProperty('')
    version = DefaultProperty('')

    from_node = DefaultProperty(Node)
    to_node = DefaultProperty(Node)

    default_connect_stack: Dict[str, str] = {}

    def __init_subclass__(cls, **kwargs):

        source = kwargs.get('source')
        cls.source.set_default(source)

        version = kwargs.get('version')
        cls.version.set_default(source)

        register = kwargs.pop('register', False)

        if register:
            from protrend.pipeline import Pipeline
            Pipeline.register_connector(cls, **kwargs)

    def __init__(self,
                 connect_stack: Dict[str, str] = None,
                 source: str = None,
                 version: str = None,
                 from_node: Type[Node] = None,
                 to_node: Type[Node] = None):
        """
        The connector object uses results obtained during the transformation procedures
        for a given neomodel node entity.
        This object is responsible for reading several files and digest/process them into relationship instances.

        A given source/database contained in the data lake must be provided
        together with the files required for the connection to take place.

        A pandas DataFrame is the main engine to read, load, process and transform data contained in these files into
        structured relationships

        :type connect_stack: Dict[str, str]
        :type source: str
        :type version: str
        :type from_node: Type[Node]
        :type to_node: Type[Node]

        :param connect_stack: Dictionary containing the pair name and file name.
        The key should be used to identify the file in the connect stack,
        whereas the value should be the file name in the data lake
        :param source: The name of the data source in the data lake (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the data lake (e.g. 0.0.0, 0.0.1, etc)
        :param from_node: The source node type associated with this connector, and thus the source of the relation.
        Note that, it should be created only one connector for each node-node relationship
        :param to_node: The target node type associated with this connector, and thus the end of the relation.
        Note that, it should be created only one connector for each node-node relationship
        """
        self._connect_stack = {}
        self._write_stack = []
        self.source = source
        self.version = version
        self.from_node = from_node
        self.to_node = to_node

        self.load_connect_stack(connect_stack)

    def load_connect_stack(self, connect_stack: Dict[str, str] = None):

        self._connect_stack = {}

        if not connect_stack:
            connect_stack = self.default_connect_stack

        for key, file in connect_stack.items():
            dl_file = os.path.join(Settings.DATA_LAKE_PATH, self.source, self.version, file)

            self._connect_stack[key] = dl_file

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def connect_stack(self) -> Dict[str, str]:
        return self._connect_stack

    @property
    def write_path(self) -> str:
        return os.path.join(Settings.DATA_LAKE_PATH, self.source, self.version)

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
    def infer_write_stack(cls) -> WriteStack:
        file_name = f'connected_{cls.from_node.default.node_name()}_{cls.to_node.default.node_name()}'
        write_stack = WriteStack(transformed=None,
                                 integrated=None,
                                 nodes=None,
                                 connected=file_name)
        return write_stack

    def stack_json(self, df: pd.DataFrame):
        name = f'connected_{self.from_node.node_name()}_{self.to_node.node_name()}'
        df = df.copy(deep=True)
        df = df.reset_index(drop=True)
        fp = os.path.join(self.write_path, f'{name}.json')
        json_partial = partial(write_json_frame, file_path=fp, df=df)
        self._write_stack.append(json_partial)

    def make_connection(self,
                        from_identifiers: List[Any],
                        to_identifiers: List[Any],
                        kwargs: dict = None) -> pd.DataFrame:

        size = len(from_identifiers)

        connection = {'from_node': [self.from_node.node_name()] * size,
                      'to_node': [self.to_node.node_name()] * size,
                      'from_identifier': from_identifiers.copy(),
                      'to_identifier': to_identifiers.copy(),
                      'load': ['create'] * size,
                      'what': ['relationships'] * size}

        if not kwargs:
            kwargs = {}

        connection.update(kwargs)

        return pd.DataFrame(connection)
