from typing import Sequence, List

from protrend.load import Loader
from protrend.log import ProtrendLogger
from protrend.transform.connector import Connector
from protrend.transform.transformer import Transformer


def sort_by_order(transformer: Transformer):
    return transformer.order


class Director:

    def __init__(self,
                 transformers: Sequence[Transformer] = None,
                 connectors: Sequence[Connector] = None,
                 loaders: Sequence[Loader] = None):

        """
        Director is responsible for managing transformers and loaders.
        It sends instructions to the transformers on how to transform a file(s) into a node or relationship
        Reading, processing, integrating and writing are executed by each transformer

        The director is also responsible for sending instructions to the loaders on how to load the transformed nodes
        into the database.

        :type transformers: Sequence[Transformer]
        :param transformers: transformers for processing staging area files into nodes

        :type connectors: Sequence[Connector]
        :param connectors: connectors for processing staging area files into relationships

        :type loaders: Sequence[Loader]
        :param loaders: loaders for converting data lake files into nodes and relationships
        """

        if not transformers:
            transformers = []

        if not connectors:
            connectors = []

        if not loaders:
            loaders = []

        self._transformers = list(transformers)
        self._connectors = list(connectors)
        self._loaders = list(loaders)

    @property
    def transformers(self) -> List[Transformer]:
        return sorted(self._transformers, key=sort_by_order, reverse=True)

    @property
    def connectors(self) -> List[Connector]:
        return self._connectors

    @property
    def loaders(self) -> List[Loader]:
        return self._loaders

    def transform(self):
        """
        Reading, processing, transforming, cleaning, integrating and writing one or more file types
        into multiple nodes.
        Transformation is performed step-wise according to the transformers order.

        :return:
        """

        transformers_info = ' '.join([transformer.node.node_name() for transformer in self.transformers])
        ProtrendLogger.log.info(f'Director call for transform in the following transformers: {transformers_info}')

        for transformer in self.transformers:

            try:

                ProtrendLogger.log.info(f'Starting transformer: {transformer}')

                df = transformer.transform()
                transformer.integrate(df)
                transformer.write()

                ProtrendLogger.log.info(f'Finishing transformer: {transformer.node.node_name()} with stats: {df.shape}')

            except:
                ProtrendLogger.log.exception(f'Exception occurred in transformer {transformer}')

    def connect(self):
        """
        Reading, processing, transforming, cleaning, integrating and writing one or more file types
        into multiple relationships.

        :return:
        """

        connectors_info = ' '.join([f'{connector.from_node.node_name()}_{connector.to_node.node_name()}\n'
                                    for connector in self.connectors])
        ProtrendLogger.log.info(f'Director call for connect in the following connectors: \n {connectors_info}')

        for connector in self.connectors:

            try:

                ProtrendLogger.log.info(f'Starting connector: {connector}')

                connector.connect()
                connector.write()

                ProtrendLogger.log.info(f'Finishing connector: '
                                        f'{connector.from_node.node_name()}_{connector.to_node.node_name()}')

            except:
                ProtrendLogger.log.exception(f'Exception occurred in connector {connector}')

    def load(self):
        """
        Loading multiple nodes and relationships. The nodes and relationships must be encoded into a csv file created
        by a transformer.
        Loading is performed step-wise according to the loaders order.

        :return:
        """

        loaders = ' '.join([f'{loader.__class__.__name__}' for loader in self.loaders])
        ProtrendLogger.log.info(f'Director call for load in the following loaders: {loaders}')

        for loader in self.loaders:
            try:
                ProtrendLogger.log.info(f'Starting loader: {loader.__class__.__name__} with the following files:')
                for file_path in loader.load_stack:
                    ProtrendLogger.log.info(f'{file_path}')

                loader.load()

                ProtrendLogger.log.info(f'Finishing loader: {loader.__class__.__name__}')

            except:
                ProtrendLogger.log.exception(f'Exception occurred in loader {loader.__class__.__name__}')
