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

    @property
    def df(self) -> pd.DataFrame:

        if self._df.empty:
            return pd.DataFrame(columns=list(self.node.cls_keys()))

        return self._df

    def read(self, *args, **kwargs):
        pass

    def process(self):
        pass

    def integrate(self, *properties):

        snapshot = self.node_snapshot()
        latest_identifier = self.node.latest_identifier()
        integer = protrend_id_decoder(latest_identifier)

        df = snapshot.query(f'name == regprecise & version == {self._version}')

        if df.empty:
            integer += 1
            regprecise = dict(protrend_id=integer,
                              name=self.name,
                              type=self.type,
                              url=self.url,
                              doi=self.doi,
                              authors=self.authors,
                              description=self.description,
                              version=self._version)

            nodes = self.node.node_from_dict(regprecise, save=True)
            node = nodes[0]
            df = node.to_df()

        self.stack_csv('source', df)
