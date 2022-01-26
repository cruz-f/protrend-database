from abc import abstractmethod

import pandas as pd

from protrend.transform import Transformer, Connector
from protrend.transform.transformations import drop_empty_string, drop_duplicates


class OperonDBTransformer(Transformer, register=False):

    @staticmethod
    def transform_operon(conserved: pd.DataFrame, known: pd.DataFrame) -> pd.DataFrame:
        conserved = conserved.dropna(subset=['coid', 'op'])
        conserved = drop_empty_string(conserved, 'coid', 'op')
        conserved = drop_duplicates(conserved, subset=['coid', 'op'], perfect_match=True)

        known = known.dropna(subset=['koid', 'op'])
        known = drop_empty_string(known, 'koid', 'op')
        known = drop_duplicates(known, subset=['koid', 'op'], perfect_match=True)

        conserved = conserved.assign(operon_db_id=conserved['coid'].copy(),
                                     operon_genes=conserved['op'].str.split(','))
        known = known.assign(operon_db_id=known['koid'].copy(),
                             operon_genes=known['op'].str.split(','))

        operon = pd.concat([conserved, known])
        operon = operon.reset_index(drop=True)
        return operon

    @abstractmethod
    def transform(self):
        pass


class OperonDBConnector(Connector, register=False):

    @abstractmethod
    def connect(self):
        pass
