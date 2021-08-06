from typing import Sequence, Set

from protrend.load.loader import Loader
from protrend.transform.transformer import Transformer


def sort_by_order(transformer: Transformer):
    return transformer.order


class Director:

    def __init__(self,
                 transformers: Sequence[Transformer] = None,
                 loaders: Sequence[Loader] = None):

        """
        Director is responsible for managing transformers and loaders.
        It sends instructions to the transformers on how to transform a file(s) into a node or relationship
        Reading, processing, integrating and writing are executed by each transformer

        The director is also responsible for sending instructions to the loaders on how to load the transformed nodes
        into the database.

        :type transformers: Sequence[Transformer]
        :param transformers: transformers for processing staging area files into nodes and relationships

        :type transformers: Sequence[Loader]
        :param transformers: loaders for converting data lake files into nodes and relationships
        """

        self._transformers = transformers
        self._loaders = loaders

    @property
    def transformers(self) -> Set[Transformer]:
        transformers = set(self._transformers)
        return set(sorted(transformers, key=sort_by_order, reverse=True))

    @property
    def loaders(self) -> Set[Loader]:
        return set(self._loaders)

    def transform(self):
        """
        Reading, processing, transforming, cleaning, integrating and writing one or more file types
        into multiple nodes and relationships.
        Transformation is performed step-wise according to the transformers order.

        :return:
        """

        for transformer in self.transformers:
            df = transformer.transform()
            transformer.integrate(df)
            transformer.write()

        for transformer in self.transformers:
            transformer.connect()
            transformer.write()

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