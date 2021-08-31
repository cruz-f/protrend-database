from typing import Sequence, List

from protrend.load.loader import Loader
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

        for transformer in self.transformers:
            df = transformer.transform()
            transformer.integrate(df)
            transformer.write()

    def connect(self):
        """
        Reading, processing, transforming, cleaning, integrating and writing one or more file types
        into multiple relationships.

        :return:
        """

        for connector in self.connectors:
            connector.connect()
            connector.write()

    def load(self):
        """
        Loading multiple nodes and relationships. The nodes and relationships must be encoded into a csv file created
        by a transformer.
        Loading is performed step-wise according to the loaders order.

        :return:
        """

        for loader in self.loaders:

            for df in loader.read():
                loader.load(df)