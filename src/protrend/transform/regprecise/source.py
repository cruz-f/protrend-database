import pandas as pd

from protrend.transform.regprecise.settings import SourceSettings
from protrend.transform.transformer import DefaultTransformer


class SourceTransformer(DefaultTransformer):
    name = 'regprecise'
    type = 'database'
    url = 'https://regprecise.lbl.gov/'
    doi = '10.1186/1471-2164-14-745'
    authors = ['Pavel S Novichkov', 'Alexey E Kazakov', 'Dmitry A Ravcheev', 'Semen A Leyn', 'Galina Y Kovaleva',
               'Roman A Sutormin', 'Marat D Kazanov', 'William Riehl', 'Adam P Arkin',
               'Inna Dubchak', 'Dmitry A Rodionov']
    description = 'RegPrecise 3.0: A resource for genome-scale exploration of transcriptional regulation in bacteria'

    default_settings = SourceSettings
    columns = {'protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'}

    def transform(self):

        regprecise = dict(name=self.name,
                          type=self.type,
                          url=self.url,
                          doi=self.doi,
                          authors=self.authors,
                          description=self.description)

        return pd.DataFrame(regprecise, index=[0])
