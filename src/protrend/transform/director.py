from typing import Tuple

from protrend.transform.transformer import Transformer


class Director:

    def __init__(self, *transformers: Transformer):

        """
        Director is responsible for managing transformers.
        It sends instructions to the transformers on how to transform a file(s) into a node or relationship
        Reading, processing, integrating and writing are executed by each transformer

        :type transformers: Transformer
        :param transformers: transformers for processing staging area files into nodes and relationships
        """

        self._transformers = transformers

    @property
    def transformers(self) -> Tuple[Transformer]:

        """
        Returns the transformers associated with this director
        """

        return self._transformers

    def transform(self):

        """
        Reading, processing, integrating and writing one or more file types into a node or relationship.
        Transformation is performed step-wise according to the transformers order.
        """

        for transformer in self._transformers:
            transformer.read()
            transformer.process()
            transformer.integrate()
            transformer.write()
