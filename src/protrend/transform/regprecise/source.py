import pandas as pd

from protrend.model.model import Source
from protrend.transform.regprecise.base import RegPreciseTransformer
from protrend.utils import SetList


class SourceTransformer(RegPreciseTransformer,
                        source='regprecise',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = 'regprecise'
    type = 'database'
    url = 'https://regprecise.lbl.gov/'
    doi = '10.1186/1471-2164-14-745'
    authors = ['Pavel S Novichkov', 'Alexey E Kazakov', 'Dmitry A Ravcheev', 'Semen A Leyn', 'Galina Y Kovaleva',
               'Roman A Sutormin', 'Marat D Kazanov', 'William Riehl', 'Adam P Arkin',
               'Inna Dubchak', 'Dmitry A Rodionov']
    description = 'RegPrecise 3.0: A resource for genome-scale exploration of transcriptional regulation in bacteria'

    columns = SetList(['name', 'type', 'url', 'doi', 'authors', 'description', 'protrend_id'])

    def transform(self):

        regprecise = dict(name=[self.name],
                          type=[self.type],
                          url=[self.url],
                          doi=[self.doi],
                          authors=[self.authors],
                          description=[self.description])

        df = pd.DataFrame(regprecise, index=[0])

        self.stack_transformed_nodes(df)

        return df
