import pandas as pd

from protrend.model import Source, Organism, Gene, Operon
from protrend.transform import SourceMixIn
from protrend.transform.operondb.base import OperonDBTransformer, OperonDBConnector
from protrend.utils import SetList


class SourceTransformer(SourceMixIn, OperonDBTransformer,
                        source='operondb',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['operondb']
    type = ['database']
    url = ['https://operondb.jp/']
    doi = ['10.1093/nar/gkq1090']
    authors = [['Shujiro Okuda', 'Akiyasu C Yoshizawa']]
    description = ['ODB: a database for operon organizations, 2011 update']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceToOrganismConnector(OperonDBConnector,
                                source='operondb',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Organism,
                                register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        df = self.create_connection(source='source', target='operon',
                                    target_column='organism', cardinality='one_to_many')
        self.stack_json(df)


class SourceConnector(OperonDBConnector,
                      source='operondb',
                      version='0.0.0',
                      register=False):

    def _connect(self, target: str) -> pd.DataFrame:
        # noinspection DuplicatedCode
        source_df, target_df = self.transform_stacks(source='source',
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={})

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        urls = []
        ids = []
        keys = []
        for _id in target_df['operon_db_id']:

            if _id.lower().startswith('k'):
                urls.append(f'https://operondb.jp/known/{_id}')
                ids.append(_id)
                keys.append('known')

            else:
                urls.append(f'https://operondb.jp/conserved/{_id}')
                ids.append(_id)
                keys.append('conserved')

        kwargs = dict(url=urls,
                      external_identifier=ids,
                      key=keys)

        return self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)


class SourceToOperonConnector(SourceConnector,
                              source='operondb',
                              version='0.0.0',
                              from_node=Source,
                              to_node=Operon,
                              register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'operon': 'integrated_operon.json'}

    def connect(self):
        df = self._connect('operon')
        self.stack_json(df)


class SourceToGeneConnector(SourceConnector,
                            source='operondb',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        df = self._connect('gene')
        self.stack_json(df)
