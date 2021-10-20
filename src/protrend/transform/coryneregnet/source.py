import pandas as pd

from protrend.model.model import Source
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.utils import SetList


class SourceTransformer(DBTBSTransformer):
    name = 'coryneregnet'
    type = 'database'
    url = 'https://www.exbio.wzw.tum.de/coryneregnet/'
    doi = '10.1038/s41597-020-0484-9'
    authors = ['Mariana Teixeira Dornelles Parise', 'Doglas Parise', 'Rodrigo Bentes Kato',
               'Josch Konstantin Pauling', 'Andreas Tauch', 'Vasco Ariston de Carvalho Azevedo', 'Jan Baumbach']
    description = 'CoryneRegNet 7, the reference database and analysis platform for corynebacterial gene regulatory networks'

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