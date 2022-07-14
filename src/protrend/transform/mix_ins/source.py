from typing import Union

import pandas as pd

from protrend.transform.transformer import Transformer
from protrend.utils import SetList


class SourceMixIn:
    name = ['']
    type = ['']
    url = ['']
    doi = ['']
    authors = [[]]
    description = ['']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])

    def transform(self: Union[Transformer, 'SourceMixIn']):
        db = dict(name=self.name,
                  type=self.type,
                  url=self.url,
                  doi=self.doi,
                  authors=self.authors,
                  description=self.description)

        df = pd.DataFrame(db, index=list(range(len(self.name))))

        self.stack_transformed_nodes(df)

        return df
