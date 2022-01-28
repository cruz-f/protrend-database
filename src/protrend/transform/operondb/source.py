from protrend.model import Source, Organism, Gene, Operon
from protrend.transform.mix_ins import SourceMixIn
from protrend.transform.operondb.base import OperonDBTransformer, OperonDBConnector
from protrend.utils import SetList
from protrend.utils.constants import DATABASE


class SourceTransformer(SourceMixIn, OperonDBTransformer,
                        source='operondb',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['operondb']
    type = [DATABASE]
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

    def connect(self):
        df = self.create_connection(source='source', target='operon',
                                    target_column='organism', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToOperonConnector(OperonDBConnector,
                              source='operondb',
                              version='0.0.0',
                              from_node=Source,
                              to_node=Operon,
                              register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target='operon',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={})

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        if 'operon_db_id' in target_df.columns:

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

            df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)

        else:
            df = self.connection_frame(source_ids=source_ids, target_ids=target_ids)

        self.stack_connections(df)


class SourceToGeneConnector(OperonDBConnector,
                            source='operondb',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):

    def connect(self):
        df = self.create_connection(source='source', target='gene',
                                    cardinality='one_to_many')
        self.stack_connections(df)
