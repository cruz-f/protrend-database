from abc import abstractmethod

import pandas as pd

from neo4j.exceptions import Neo4jError, DriverError

from protrend.transform import Transformer


class StandardizerTransformer(Transformer, register=False):

    @abstractmethod
    def transform(self):
        pass

    def fetch_nodes(self):
        try:
            return self.node.nodes.all()

        except (Neo4jError, DriverError):
            return []

    def integrate(self, df: pd.DataFrame):
        df.assign(load='delete', what='nodes')
        self.stack_integrated_nodes(df)
        self.stack_nodes(df)
