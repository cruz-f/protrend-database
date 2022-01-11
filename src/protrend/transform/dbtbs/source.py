from protrend.transform import SourceMixIn
from protrend.model import (Source, Organism, RegulatoryFamily, Regulator, Gene, TFBS,
                            RegulatoryInteraction)
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
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
    default_connect_stack = {'source': 'integrated_source.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        df = self.create_connection(source='source', target='organism', cardinality='one_to_many')
        self.stack_json(df)


class SourceToRegulatoryFamilyConnector(DBTBSConnector,
                                        source='dbtbs',
                                        version='0.0.4',
                                        from_node=Source,
                                        to_node=RegulatoryFamily,
                                        register=True):
    default_connect_stack = {'source': 'integrated_source.json',
                             'rfam': 'integrated_regulatoryfamily.json'}

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
        self.stack_json(df)


class SourceConnector(DBTBSConnector, source='dbtbs', version='0.0.4', register=False):

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
    default_connect_stack = {'source': 'integrated_source.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        df = self._connect(target='regulator', external_id_col='name_dbtbs', external_url_col='url', key_id='tfac')
        self.stack_json(df)


class SourceToGeneConnector(SourceConnector,
                            source='dbtbs',
                            version='0.0.4',
                            from_node=Source,
                            to_node=Gene,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        df = self._connect(target='gene', external_id_col='tf', external_url_col='url', key_id='tfac')
        self.stack_json(df)


class SourceToTFBSConnector(SourceConnector,
                            source='dbtbs',
                            version='0.0.4',
                            from_node=Source,
                            to_node=TFBS,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self._connect(target='tfbs', external_id_col='tf', external_url_col='url', key_id='tfac')
        self.stack_json(df)


class SourceToRegulatoryInteractionConnector(SourceConnector,
                                             source='dbtbs',
                                             version='0.0.4',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):
    default_connect_stack = {'source': 'integrated_source.json',
                             'ri': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self._connect(target='ri', external_id_col='tf', external_url_col='url', key_id='tfac')
        self.stack_json(df)
