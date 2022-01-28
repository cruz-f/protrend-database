import os
from typing import List

import pandas as pd
from neomodel import AttemptedCardinalityViolation
from tqdm import tqdm

from protrend.io import read_json_frame
from protrend.log import ProtrendLogger
from protrend.model import get_node_by_name, get_nodes_relationships, connect_nodes
from protrend.utils import DefaultProperty, build_file_path


def load_stack_from_source_version(source, version) -> List[str]:
    load_stack = []

    from protrend.pipeline import Pipeline
    transformers = Pipeline.default_transformers.get((source, version), [])
    connectors = Pipeline.default_connectors.get((source, version), [])

    for transformer in transformers:
        file = transformer.nodes_file()
        nodes = build_file_path(source=source, version=version, file=file)
        load_stack.append(nodes)

    for connector in connectors:
        file = connector.connected_file()
        connections = build_file_path(source=source, version=version, file=file)
        load_stack.append(connections)

    return load_stack


class Loader:
    source = DefaultProperty()
    version = DefaultProperty()

    def __init_subclass__(cls, **kwargs):

        source = kwargs.get('source')
        cls.source.set_default(cls, source)

        version = kwargs.get('version')
        cls.version.set_default(cls, version)

        register = kwargs.pop('register', False)

        if register:
            from protrend.pipeline import Pipeline
            Pipeline.register_loader(cls, **kwargs)

    def __init__(self,
                 source: str = None,
                 version: str = None):

        """
        The load object must implements loading procedures for a set of data lake files.
        Using these files, the loader will create or update existing nodes in the neo4j database.

        :type source: str
        :type version: str

        :param source: The name of the data source in the data lake (e.g. regprecise, collectf, etc)
        :param version: The version of the data source in the data lake (e.g. 0.0.0, 0.0.1, etc)
        """
        self.source = source
        self.version = version

    # --------------------------------------------------------
    # Loader API
    # --------------------------------------------------------
    def load(self):

        """
        The method responsible for reading files into pandas DataFrames
        and creating or updating nodes and relationships in the database

        :return:
        """

        load_stack = load_stack_from_source_version(self.source, self.version)

        for file_path in load_stack:

            if os.path.exists(file_path):
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

        from_node_name = df['from_node'].iloc[0]
        from_node = get_node_by_name(from_node_name)
        from_nodes = from_node.node_to_dict(to='node')

        to_node_name = df['to_node'].iloc[0]
        to_node = get_node_by_name(to_node_name)
        to_nodes = to_node.node_to_dict(to='node')

        from_rels, to_rels = get_nodes_relationships(from_node=from_node, to_node=to_node)

        for _, relationship in tqdm(df.iterrows(),
                                    desc=f'relationship: {from_node_name} - {to_node_name}',
                                    total=df.shape[0]):

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
