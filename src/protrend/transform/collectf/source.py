from protrend.model import (Source, Organism, RegulatoryFamily, Regulator, Gene, TFBS, RegulatoryInteraction)
from protrend.transform.collectf.base import CollecTFTransformer, CollecTFConnector
from protrend.transform.mix_ins import SourceMixIn
from protrend.utils import SetList, is_null
from protrend.utils.processors import to_list_nan


class SourceTransformer(SourceMixIn, CollecTFTransformer,
                        source='collectf',
                        version='0.0.1',
                        node=Source,
                        order=100,
                        register=True):
    name = ['collectf']
    type = ['database']
    url = ['http://collectf.org/']
    doi = ['10.1093/nar/gkt1123']
    authors = [['Sefa Kiliç', 'Elliot R White', 'Dinara M Sagitova', 'Joseph P Cornish', 'Ivan Erill']]
    description = ['CollecTF: a database of experimentally validated transcription factor-binding sites in Bacteria']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceToOrganismConnector(CollecTFConnector,
                                source='collectf',
                                version='0.0.1',
                                from_node=Source,
                                to_node=Organism,
                                register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target='organism',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={})

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        url = []
        ext_id = []
        key = []
        if 'taxonomy' in target_df.columns:
            for tax_id in target_df['taxonomy']:
                if is_null(tax_id):
                    url.append(None)
                    ext_id.append(None)
                    key.append(None)
                else:
                    url.append(f'http://www.collectf.org/browse/view_motif_reports_by_taxonomy/{tax_id}')
                    ext_id.append(tax_id)
                    key.append('view_motif_reports_by_taxonomy')

        kwargs = dict(url=url,
                      external_identifier=ext_id,
                      key=key)

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)
        self.stack_connections(df)


class SourceToRegulatoryFamilyConnector(CollecTFConnector,
                                        source='collectf',
                                        version='0.0.1',
                                        from_node=Source,
                                        to_node=RegulatoryFamily,
                                        register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'rfam': 'integrated_regulatoryfamily.json'}

    def connect(self):
        df = self.create_connection(source='source',
                                    target='rfam',
                                    cardinality='one_to_many')
        self.stack_connections(df)


class SourceToRegulatorConnector(CollecTFConnector,
                                 source='collectf',
                                 version='0.0.1',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target='regulator',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={})

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        url = []
        ext_id = []
        key = []
        if 'uniprot_accession' in target_df:
            for reg_id in target_df['uniprot_accession']:
                if is_null(reg_id):
                    url.append(None)
                    ext_id.append(None)
                    key.append(None)
                else:
                    url.append(f'http://www.collectf.org/uniprot/{reg_id}')
                    ext_id.append(reg_id)
                    key.append('uniprot')

        kwargs = dict(url=url,
                      external_identifier=ext_id,
                      key=key)

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)
        self.stack_connections(df)


class SourceToGeneConnector(CollecTFConnector,
                            source='collectf',
                            version='0.0.1',
                            from_node=Source,
                            to_node=Gene,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        df = self.create_connection(source='source',
                                    target='gene',
                                    cardinality='one_to_many')
        self.stack_connections(df)


class SourceToTFBSConnector(CollecTFConnector,
                            source='collectf',
                            version='0.0.1',
                            from_node=Source,
                            to_node=TFBS,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self.create_connection(source='source',
                                    target='tfbs',
                                    cardinality='one_to_many')
        self.stack_connections(df)


class SourceToRegulatoryInteractionConnector(CollecTFConnector,
                                             source='collectf',
                                             version='0.0.1',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target='rin',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={'regulon': [to_list_nan]})

        target_df = target_df.explode('regulon')

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        url = []
        ext_id = []
        key = []
        for reg_id in target_df['regulon']:
            if is_null(reg_id):
                url.append(None)
                ext_id.append(None)
                key.append(None)
            else:
                url.append(f'http://www.collectf.org/uniprot/{reg_id}')
                ext_id.append(reg_id)
                key.append('uniprot')

        kwargs = dict(url=url,
                      external_identifier=ext_id,
                      key=key)

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)
        self.stack_connections(df)
