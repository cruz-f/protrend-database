import os
from abc import ABCMeta, abstractmethod
from functools import partial
from typing import List, Any, Type, Callable, Dict, Set

import pandas as pd

from protrend.model.node import Node
from protrend.transform.settings import ConnectorSettings
from protrend.utils.settings import DATA_LAKE_PATH


class Connector(metaclass=ABCMeta):
    """
    Connector interface.

    The following methods must be implemented to set up a connector for each relationship

    A connector is responsible for reading, processing, integrating and writing relationships files.

    A connector starts with data from the staging area and ends with structured relationships.
    """

    def __init__(self, settings: ConnectorSettings):
        """
        The connector object contains and uses transformation procedures for a given neomodel node entity.
        This object is responsible for reading several files and digest/process them into relationship instances.

        A given source/database contained in the staging area must be provided
        together with the files required for the transformation to take place.

        A pandas DataFrame is the main engine to read, load, process and transform data contained in these files into
        structured relationships

        :param settings: a ConnectorSettings than contains all settings for source,
        version and files to perform the connection on
        """

        self._settings = settings

        self._connect_stack = {}
        self._write_stack = []

        self._load_settings()

    def _load_settings(self):
        for key, file in self._settings.connect.items():
            file_path = os.path.join(DATA_LAKE_PATH, self.source, self.version, file)

            self._connect_stack[key] = file_path

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def settings(self) -> ConnectorSettings:
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
    def from_node(self) -> Type[Node]:
        return self.settings.from_node

    @property
    def to_node(self) -> Type[Node]:
        return self.settings.to_node

    @property
    def connect_stack(self) -> Dict[str, str]:
        return self._connect_stack

    @property
    def write_stack(self) -> List[Callable]:
        return self._write_stack

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

    def __eq__(self, other: 'Connector'):
        return self.__class__.__name__ and other.__class__.__name__ and self.settings == other.settings

    @abstractmethod
    def connect(self):
        """
        The method responsible for connecting an integrated pandas DataFrames to the remaining database.
        The connect method should set up multiple connections into several pandas DataFrames
        using either protrend identifiers or node factors

        Interface implementation with unknown signature.
        Concrete implementations are available at the transformer children.

        :return:
        """
        pass

    def write(self):

        if not os.path.exists(self.write_path):
            os.makedirs(self.write_path)

        for csv in self._write_stack:
            csv()

        self._write_stack = []

    def stack_csv(self, df: pd.DataFrame):

        name = f'connected_{self.from_node.node_name()}_{self.to_node.node_name()}'

        df_copy = df.copy(deep=True)
        fp = os.path.join(self.write_path, f'{name}.csv')
        csv = partial(df_copy.to_csv, path_or_buf=fp)
        self._write_stack.append(csv)

    def make_connection(self,
                        size: int,
                        from_identifiers: List[Any],
                        to_identifiers: List[Any],
                        kwargs: dict = None) -> pd.DataFrame:

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


class DefaultConnector(Connector):
    default_settings: Type[ConnectorSettings] = ConnectorSettings

    def __init__(self, settings: ConnectorSettings = None):
        if not settings:
            settings = self.default_settings()

        super().__init__(settings)

    @abstractmethod
    def connect(self):
        pass
