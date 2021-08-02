import os
from abc import abstractmethod
from typing import Dict, Tuple

import pandas as pd

from protrend.load.settings import LoaderSettings
from protrend.model.node import Node
from protrend.transform.settings import TransformerSettings
from protrend.utils.settings import DATA_LAKE_PATH


class Loader:

    def __init__(self, settings: LoaderSettings):

        """

        :param settings: a TransformerSettings than contains all settings for source,
        version and files to perform the transformation on
        """

        self._settings = settings
        self._files = {}
        self._attrs = {}

        self._load_settings()

    def _load_settings(self):

        for key, file in self._settings.files.items():

            file_path = os.path.join(DATA_LAKE_PATH, self.source, self.version, file)

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
    def read_path(self) -> str:
        return os.path.join(DATA_LAKE_PATH, self.source, self.version)

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

    def load(self, df: pd.DataFrame):

        # update nodes
        self._nodes_to_update(df=df)

        # create nodes
        self._nodes_to_create(df=df)

        # create nodes
        self._relationships_to_create(df=df)

    def _nodes_to_update(self, df: pd.DataFrame) -> pd.DataFrame:

        mask = (df['load'] == 'update') & (df['what'] == 'nodes')
        df = df.loc[mask, :]

        if df.empty:
            return

        self.node.node_from_df(nodes=df, save=True)

    def _nodes_to_create(self, df: pd.DataFrame) -> pd.DataFrame:

        mask = (df['load'] == 'create') & (df['what'] == 'nodes')
        df = df.loc[mask, :]

        if df.empty:
            return

        self.node.node_update_from_df(nodes=df, save=True)

    def _relationships_to_create(self, df: pd.DataFrame) -> pd.DataFrame:

        mask = (df['load'] == 'create') & (df['what'] == 'relationships')
        df = df.loc[mask, :]

        if df.empty:
            return

        self.node.node_update_from_df(nodes=df, save=True)
