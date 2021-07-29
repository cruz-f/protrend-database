from typing import Tuple, Dict

import pandas as pd

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
        self._relationships = {}

    @property
    def transformers(self) -> Tuple[Transformer]:

        """
        Returns the transformers associated with this director
        """

        return self._transformers

    @property
    def relationships(self) -> Dict[str, pd.DataFrame]:

        """
        Returns the relationships associated with this director
        """

        return self._relationships

    def connect(self, from_node: str, to_node: str, df: pd.DataFrame):
        pass

    def integrate(self):

        """
        Reading, processing, integrating and writing one or more file types into a node or relationship.
        Transformation is performed step-wise according to the transformers order.
        """

        for transformer in self._transformers:
            transformer.read()
            df = transformer.transform()
            df = transformer.load_nodes(df)
            relationships_dict = transformer.load_relationships(df)
            self.relationships.update(relationships_dict)
            transformer.write()

        for (from_node, to_node), df in self.relationships.items():
            self.connect(from_node, to_node, df)
