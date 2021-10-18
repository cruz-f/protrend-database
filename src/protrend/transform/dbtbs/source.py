import pandas as pd

from protrend.model.model import Source
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.utils import SetList


class SourceTransformer(DBTBSTransformer):
    name = 'dbtbs'
    type = 'database'
    url = 'https://dbtbs.hgc.jp/'
    doi = '10.1093/nar/gkm910'
    authors = ['Nicolas Sierro', 'Yuko Makita', 'Michiel de Hoon', 'Kenta Nakai']
    description = 'DBTBS: a database of transcriptional regulation in Bacillus subtilis containing upstream intergenic conservation information'

    default_node = Source
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])

    def transform(self):
        db = dict(name=[self.name],
                        type=[self.type],
                        url=[self.url],
                        doi=[self.doi],
                        authors=[self.authors],
                        description=[self.description])

        df = pd.DataFrame(db, index=[0])

        self._stack_transformed_nodes(df)

        return df