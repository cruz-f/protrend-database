import os
from typing import List

import pandas as pd

from protrend.io.csv import read_csv
from protrend.load.settings import LoaderSettings
from protrend.model.node import Node
from protrend.utils.settings import DATA_LAKE_PATH


class Loader:

    def __init__(self, settings: LoaderSettings):

        """

        :param settings: a TransformerSettings than contains all settings for source,
        version and files to perform the transformation on
        """

        self._settings = settings
        self._files = []

        self._load_settings()

    def _load_settings(self):

        for file in self._settings.files:

            file_path = os.path.join(DATA_LAKE_PATH, self.source, self.version, file)

            if os.path.exists(file_path):

                self._files.append(file_path)
            else:
                raise FileNotFoundError(f'Could not found file {file_path}')

    # --------------------------------------------------------
    # Static properties
    # --------------------------------------------------------
    @property
    def settings(self) -> LoaderSettings:
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
    def files(self) -> List[str]:
        return self._files

    @property
    def node(self) -> Node:
        return self.settings.node

    @property
    def read_path(self) -> str:
        return os.path.join(DATA_LAKE_PATH, self.source, self.version)

    # --------------------------------------------------------
    # Transformer API
    # --------------------------------------------------------
    def read(self, **kwargs):

        """
        The method responsible for reading files into pandas DataFrames.
        The pandas DataFrames are yield as result.

        :param kwargs: kwargs for the pandas read_csv API
        :return: generator of pandas DataFrame
        """

        for file_path in self._files:
            yield read_csv(file_path, **kwargs)

    def load(self, df: pd.DataFrame):

        # update nodes
        self._nodes_to_update(df=df)

        # create nodes
        self._nodes_to_create(df=df)

        # create nodes
        self._relationships_to_create(df=df)

    def _nodes_to_update(self, df: pd.DataFrame):

        mask = (df['load'] == 'update') & (df['what'] == 'nodes')
        df = df.loc[mask, :]

        if df.empty:
            return

        self.node.node_from_df(nodes=df, save=True)

    def _nodes_to_create(self, df: pd.DataFrame):

        mask = (df['load'] == 'create') & (df['what'] == 'nodes')
        df = df.loc[mask, :]

        if df.empty:
            return

        self.node.node_update_from_df(nodes=df, save=True)

    def _relationships_to_create(self, df: pd.DataFrame):

        mask = (df['load'] == 'create') & (df['what'] == 'relationships')
        df = df.loc[mask, :]

        if df.empty:
            return

        self.node.node_update_from_df(nodes=df, save=True)
