from typing import Dict

import pandas as pd

from protrend.model.model import Source
from protrend.model.node import protrend_id_decoder
from protrend.transform.regprecise.settings import RegPreciseTransformSettings
from protrend.transform.transformer import Transformer


class SourceTransformer(Transformer):
    node = Source
    name = 'regprecise'
    type = 'database'
    url = ''
    doi = ''
    authors = []
    description = ''

    def __init__(self,
                 source: str = None,
                 version: str = None,
                 **files: Dict[str, str]):

        if not source:
            source = RegPreciseTransformSettings.source

        if not version:
            version = RegPreciseTransformSettings.version

        super().__init__(source=source, version=version, **files)

    def read(self, *args, **kwargs):
        pass

    def transform(self):
        pass

    def load(self, *properties):

        snapshot = self.node_snapshot()

        last_node = self.node.last_node()
        if last_node is None:
            integer = 0

        else:
            integer = protrend_id_decoder(last_node.protrend_id)

        df = snapshot.query(f'name == {self.name}')

        if df.empty:
            integer += 1
            regprecise = dict(protrend_id=integer,
                              name=self.name,
                              type=self.type,
                              url=self.url,
                              doi=self.doi,
                              authors=self.authors,
                              description=self.description)

            self.node.node_from_dict(regprecise, save=True)

            df = pd.DataFrame(regprecise, index=[0])

        self.stack_csv('source', df)
