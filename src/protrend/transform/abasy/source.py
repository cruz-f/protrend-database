import pandas as pd

from protrend.model.model import Source
from protrend.transform.abasy.base import AbasyTransformer
from protrend.utils import SetList


class SourceTransformer(AbasyTransformer):
    name = 'abasy'
    type = 'database'
    url = 'https://abasy.ccg.unam.mx'
    doi = '10.1016/j.csbj.2020.05.015'
    authors = ['Juan M. Escorcia-Rodríguez, AndreasTauch, Julio A. Freyre-Gonzáleza']
    description = 'Abasy Atlas v2.2: The most comprehensive and up-to-date inventory of meta-curated, historical, bacterial regulatory networks, their completeness and system-level characterization'

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
