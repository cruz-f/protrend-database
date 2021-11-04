import os
from typing import List

import pandas as pd
from neomodel import AttemptedCardinalityViolation

from protrend.io.json import read_json_frame
from protrend.log import ProtrendLogger
from protrend.model.node import get_node_by_name, get_nodes_relationships, connect_nodes
from protrend.utils import Settings


class Loader:
    default_load_stack: List[str] = []
    default_source: str = ''
    default_version: str = '0.0.0'

    def __init__(self,
                 load_stack: List[str] = None,
                 source: str = None,
                 version: str = None):

        """
        The load object must implements loading procedures for a set of data lake files.
        Using these files, the loader will create or update existing nodes in the neo4j database.

        :type load_stack: List[str]
        :type source: str
        :type version: str

        :param load_stack: List containing the files. The value should be the file name in the data lake
        :param source: The name of the data source in the data lake (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the data lake (e.g. 0.0.0, 0.0.1, etc)
        """
        self._load_stack = []
        self._source = source
        self._version = version

        self.load_transform_stack(load_stack)

    def load_transform_stack(self, load_stack: List[str] = None):

        self._load_stack = []

        if not load_stack:
            load_stack = self.default_load_stack

        for file in load_stack:

            dl_file = os.path.join(Settings.DATA_LAKE_PATH, self.source, self.version, file)

            if os.path.exists(dl_file):
                self._load_stack.append(dl_file)

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
    def load_stack(self) -> List[str]:
        return self._load_stack

    @property
    def read_path(self) -> str:
        return os.path.join(Settings.DATA_LAKE_PATH, self.source, self.version)

    # --------------------------------------------------------
    # Transformer API
    # --------------------------------------------------------
    def load(self):

        """
        The method responsible for reading files into pandas DataFrames
        and creating or updating nodes and relationships in the database

        :return:
        """

        for file_path in self.load_stack:
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

                except AttemptedCardinalityViolation:
                    continue

                except:
                    ProtrendLogger.log.exception(f'Could not connect {from_node_instance.identifier} to '
                                                 f'{to_node_instance.identifier} using {relationship} relation')

            # reverse connection
            for attr in to_rels:

                try:
                    connect_nodes(from_node=to_node_instance, to_node=from_node_instance, relationship=attr,
                                  kwargs=relationship)

                except AttemptedCardinalityViolation:
                    continue

                except:
                    ProtrendLogger.log.exception(f'Could not connect {to_node_instance.identifier} to '
                                                 f'{from_node_instance.identifier} using {relationship} relation')
