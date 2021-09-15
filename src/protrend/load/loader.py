import os
from typing import List, Type

import pandas as pd
from neomodel import RelationshipManager

from protrend.io.json import read_json_frame
from protrend.load.settings import LoaderSettings
from protrend.log.logger import ProtrendLogger
from protrend.model.node import get_node_by_name, get_nodes_relationships, connect_nodes
from protrend.utils.settings import DATA_LAKE_PATH


class Loader:
    default_settings: Type[LoaderSettings] = LoaderSettings

    def __init__(self, settings: LoaderSettings = None):

        """

        :param settings: a TransformerSettings than contains all settings for source,
        version and files to perform the transformation on
        """
        if not settings:
            settings = self.default_settings()

        self._settings = settings
        self._files = []

        self._load_settings()

    def _load_settings(self):

        for file in self._settings.files:

            file_path = os.path.join(DATA_LAKE_PATH, self.source, self.version, file)

            if os.path.exists(file_path):
                self._files.append(file_path)

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
    def read_path(self) -> str:
        return os.path.join(DATA_LAKE_PATH, self.source, self.version)

    # --------------------------------------------------------
    # Transformer API
    # --------------------------------------------------------
    def load(self):

        """
        The method responsible for reading files into pandas DataFrames
        and creating or updating nodes and relationships in the database

        :return:
        """

        for file_path in self._files:
            df = read_json_frame(file_path)

            # update nodes
            self._nodes_to_update(df=df)

            # create nodes
            self._nodes_to_create(df=df)

            # create nodes
            self._relationships_to_create(df=df)

    @staticmethod
    def _nodes_to_update(df: pd.DataFrame):

        mask = (df['load'] == 'update') & (df['what'] == 'nodes')
        df = df.loc[mask, :]

        if df.empty:
            return

        node_name = df['node'].iloc[0]
        node = get_node_by_name(node_name)

        node.node_update_from_df(nodes=df, save=True)

    @staticmethod
    def _nodes_to_create(df: pd.DataFrame):

        mask = (df['load'] == 'create') & (df['what'] == 'nodes')
        df = df.loc[mask, :]

        if df.empty:
            return

        node_name = df['node'].iloc[0]
        node = get_node_by_name(node_name)

        node.node_from_df(nodes=df, save=True)

    @staticmethod
    def _relationships_to_create(df: pd.DataFrame):

        mask = (df['load'] == 'create') & (df['what'] == 'relationships')
        df = df.loc[mask, :]

        if df.empty:
            return

        from_node = df['from_node'].iloc[0]
        from_node = get_node_by_name(from_node)
        from_nodes = from_node.node_to_dict(to='node')

        to_node = df['to_node'].iloc[0]
        to_node = get_node_by_name(to_node)
        to_nodes = to_node.node_to_dict(to='node')

        from_rels, to_rels = get_nodes_relationships(from_node=from_node, to_node=to_node)

        for _, relationship in df.iterrows():

            relationship = relationship.to_dict()

            from_id = relationship['from_identifier']
            from_node_instance = from_nodes.get(from_id)
            to_id = relationship['to_identifier']
            to_node_instance = to_nodes.get(to_id)

            if not from_node_instance or not to_node_instance:
                continue

            # forward connection
            for attr in from_rels:

                try:
                    connect_nodes(from_node=from_node_instance, to_node=to_node_instance, relationship=attr,
                                  kwargs=relationship)

                except:
                    ProtrendLogger.log.exception(f'Could not connect {from_node_instance.identifier} to '
                                         f'{to_node_instance.identifier} using {relationship} relation')

            # reverse connection
            for attr in to_rels:

                try:
                    connect_nodes(from_node=to_node_instance, to_node=from_node_instance, relationship=attr,
                                  kwargs=relationship)

                except:
                    ProtrendLogger.log.exception(f'Could not connect {to_node_instance.identifier} to '
                                         f'{from_node_instance.identifier} using {relationship} relation')
