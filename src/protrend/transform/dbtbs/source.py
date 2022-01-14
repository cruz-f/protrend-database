from protrend.model import (Source, Organism, RegulatoryFamily, Regulator, Gene, TFBS,
                            RegulatoryInteraction)
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.mix_ins import SourceMixIn
from protrend.utils import SetList, is_null


class SourceTransformer(SourceMixIn, DBTBSTransformer,
                        source='dbtbs',
                        version='0.0.4',
                        node=Source,
                        order=100,
                        register=True):
    name = ['dbtbs']
    type = ['database']
    url = ['https://dbtbs.hgc.jp/']
    doi = ['10.1093/nar/gkm910']
    authors = [['Nicolas Sierro', 'Yuko Makita', 'Michiel de Hoon', 'Kenta Nakai']]
    description = ['DBTBS: a database of transcriptional regulation in Bacillus subtilis containing upstream '
                   'intergenic conservation information']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceToOrganismConnector(DBTBSConnector,
                                source='dbtbs',
                                version='0.0.4',
                                from_node=Source,
                                to_node=Organism,
                                register=True):

    def connect(self):
        df = self.create_connection(source='source', target='organism', cardinality='one_to_many')
        self.stack_connections(df)


class SourceToRegulatoryFamilyConnector(DBTBSConnector,
                                        source='dbtbs',
                                        version='0.0.4',
                                        from_node=Source,
                                        to_node=RegulatoryFamily,
                                        register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target='rfam',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={})

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        url = []
        ext_id = []
        key = []
        if 'name' in target_df.columns:
            for rfam_name in target_df['name']:
                if not is_null(rfam_name) and rfam_name == 'Helix turn helix family':
                    url.append('https://dbtbs.hgc.jp/tfactable.html#HTH')
                    ext_id.append('HTH')
                    key.append('tfactable.html#')
                elif not is_null(rfam_name) and rfam_name == 'Sigma factors':
                    url.append('https://dbtbs.hgc.jp/tfactable.html#sigma')
                    ext_id.append('sigma')
                    key.append('tfactable.html#')
                else:
                    url.append('https://dbtbs.hgc.jp/')
                    ext_id.append(None)
                    key.append(None)

        kwargs = dict(url=url,
                      external_identifier=ext_id,
                      key=key)

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)
        self.stack_connections(df)


class SourceConnector(DBTBSConnector, register=False):

    def _connect(self, target: str, external_id_col: str, external_url_col: str, key_id: str):
        # noinspection DuplicatedCode
        source_df, target_df = self.transform_stacks(source='source',
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={})

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        url = []
        ext_id = []
        key = []
        if external_id_col in target_df.columns and external_url_col in target_df.columns:
            for reg_id, reg_url in zip(target_df[external_id_col], target_df[external_url_col]):
                if not is_null(reg_id) and not is_null(reg_url):
                    url.append(reg_url)
                    ext_id.append(reg_id)
                    key.append(key_id)
                else:
                    url.append(None)
                    ext_id.append(None)
                    key.append(None)

        kwargs = dict(url=url,
                      external_identifier=ext_id,
                      key=key)

        return self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)


class SourceToRegulatorConnector(SourceConnector,
                                 source='dbtbs',
                                 version='0.0.4',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):

    def connect(self):
        df = self._connect(target='regulator', external_id_col='dbtbs_name', external_url_col='url', key_id='tfac')
        self.stack_connections(df)


class SourceToGeneConnector(SourceConnector,
                            source='dbtbs',
                            version='0.0.4',
                            from_node=Source,
                            to_node=Gene,
                            register=True):

    def connect(self):
        df = self._connect(target='gene', external_id_col='tf', external_url_col='url', key_id='tfac')
        self.stack_connections(df)


class SourceToTFBSConnector(SourceConnector,
                            source='dbtbs',
                            version='0.0.4',
                            from_node=Source,
                            to_node=TFBS,
                            register=True):

    def connect(self):
        df = self._connect(target='tfbs', external_id_col='tf', external_url_col='url', key_id='tfac')
        self.stack_connections(df)


class SourceToRegulatoryInteractionConnector(SourceConnector,
                                             source='dbtbs',
                                             version='0.0.4',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):

    def connect(self):
        df = self._connect(target='rin', external_id_col='tf', external_url_col='url', key_id='tfac')
        self.stack_connections(df)
