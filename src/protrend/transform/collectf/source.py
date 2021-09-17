import pandas as pd

from protrend.model.model import Source
from protrend.transform.collectf.base import CollectfTransformer


class SourceTransformer(CollectfTransformer):
    name = 'collectf'
    type = 'database'
    url = 'http://collectf.org/'
    doi = '10.1093/nar/gkt1123'
    authors = ['Sefa Kiliç', 'Elliot R White', 'Dinara M Sagitova', 'Joseph P Cornish', 'Ivan Erill']
    description = 'CollecTF: a database of experimentally validated transcription factor-binding sites in Bacteria'

    default_node = Source
    default_node_factors = ('name',)
    default_order = 100
    columns = {'protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'}

    def transform(self):
        collectf = dict(name=[self.name],
                        type=[self.type],
                        url=[self.url],
                        doi=[self.doi],
                        authors=[self.authors],
                        description=[self.description])

        df = pd.DataFrame(collectf, index=[0])

        self._stack_transformed_nodes(df)

        return df
