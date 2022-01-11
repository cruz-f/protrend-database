from typing import Dict, Callable, List

from protrend.io import read_from_stack, read_json_frame
from protrend.model.model import (Source, Effector, Gene, Organism, Pathway, Regulator, RegulatoryFamily,
                                  RegulatoryInteraction, TFBS)
from protrend.transform import SourceMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.regulatory_family import RegulatoryFamilyTransformer
from protrend.utils import SetList
from protrend.utils.processors import to_list_nan


class SourceTransformer(SourceMixIn, RegPreciseTransformer,
                        source='regprecise',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['regprecise']
    type = ['database']
    url = ['https://regprecise.lbl.gov/']
    doi = ['10.1186/1471-2164-14-745']
    authors = [['Pavel S Novichkov', 'Alexey E Kazakov', 'Dmitry A Ravcheev', 'Semen A Leyn', 'Galina Y Kovaleva',
                'Roman A Sutormin', 'Marat D Kazanov', 'William Riehl', 'Adam P Arkin',
                'Inna Dubchak', 'Dmitry A Rodionov']]
    description = ['RegPrecise 3.0: A resource for genome-scale exploration of transcriptional regulation in bacteria']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceConnector(RegPreciseConnector,
                      source='regprecise',
                      version='0.0.0',
                      register=False):

    def _connect(self,
                 target: str,
                 target_processors: Dict[str, List[Callable]],
                 url: str,
                 external_identifier: str,
                 key: str,
                 explode: str = None):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors=target_processors)

        if explode:
            target_df = target_df.explode(explode)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        size = len(target_ids)
        kwargs = dict(url=target_df[url].to_list(),
                      external_identifier=target_df[external_identifier].to_list(),
                      key=[key] * size)

        return self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)


class SourceToEffectorConnector(SourceConnector,
                                source='regprecise',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Effector,
                                register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'effector': 'integrated_effector.json'}

    def connect(self):
        df = self._connect(target='effector', target_processors={},
                           url='url', external_identifier='effector_id', key='effector_id')
        self.stack_json(df)


class SourceToGeneConnector(SourceConnector,
                            source='regprecise',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'gene': 'integrated_gene.json'}

    def connect(self):
        df = self._connect(target='gene', target_processors={'regulon': [to_list_nan]},
                           url='url', external_identifier='regulon', key='effector_id',
                           explode='regulon')
        self.stack_json(df)


class SourceToOrganismConnector(SourceConnector,
                                source='regprecise',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Organism,
                                register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'organism': 'integrated_organism.json'}

    def connect(self):
        df = self._connect(target='organism', target_processors={},
                           url='url', external_identifier='genome_id', key='genome_id')
        self.stack_json(df)


class SourceToPathwayConnector(SourceConnector,
                               source='regprecise',
                               version='0.0.0',
                               from_node=Source,
                               to_node=Pathway,
                               register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'pathway': 'integrated_pathway.json'}

    def connect(self):
        df = self._connect(target='pathway', target_processors={},
                           url='url', external_identifier='pathway_id', key='pathway_id')
        self.stack_json(df)


class SourceToRegulatorConnector(SourceConnector,
                                 source='regprecise',
                                 version='0.0.0',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'regulator': 'integrated_regulator.json'}

    def connect(self):
        df = self._connect(target='regulator', target_processors={},
                           url='url', external_identifier='regulon_id', key='regulon_id')
        self.stack_json(df)


class SourceToRegulatoryFamilyConnector(RegPreciseConnector,
                                        source='regprecise',
                                        version='0.0.0',
                                        from_node=Source,
                                        to_node=RegulatoryFamily,
                                        register=True):
    default_connect_stack = {'source': 'integrated_source.json',
                             'rfam': 'integrated_regulatoryfamily.json'}

    def connect(self):
        source = read_from_stack(stack=self.connect_stack, key='source',
                                 columns=SourceTransformer.columns, reader=read_json_frame)
        target = read_from_stack(stack=self.connect_stack, key='rfam',
                                 columns=RegulatoryFamilyTransformer.columns, reader=read_json_frame)

        protrend_id = source.loc[0, 'protrend_id']

        from_identifiers = []
        to_identifiers = []
        urls = []
        external_ids = []
        keys = []
        for target_key in ('tffamily_id', 'collection_id', 'riboswitch_id'):
            target_df = target.dropna(subset=[target_key])

            from_identifier = target_df['protrend_id'].to_list()

            size = len(from_identifier)
            to_identifier = [protrend_id] * size

            url = target_df['url'].to_list()
            external_id = target_df[target_key].to_list()
            keyword = [target_key] * size

            from_identifiers.extend(from_identifier)
            to_identifiers.extend(to_identifier)
            urls.extend(url)
            external_ids.extend(external_id)
            keys.extend(keyword)

        kwargs = dict(url=urls, external_identifier=external_ids, key=keys)

        df = self.connection_frame(source_ids=from_identifiers, target_ids=to_identifiers, kwargs=kwargs)
        self.stack_json(df)


class SourceToRegulatoryInteractionConnector(SourceConnector,
                                             source='regprecise',
                                             version='0.0.0',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self._connect(target='rin', target_processors={},
                           url='url', external_identifier='regulon_id', key='regulon_id')
        self.stack_json(df)


class SourceToTFBSConnector(SourceConnector,
                            source='regprecise',
                            version='0.0.0',
                            from_node=Source,
                            to_node=TFBS,
                            register=True):
    default_connect_stack = {'source': 'integrated_source.json', 'tfbs': 'integrated_tfbs.json'}

    def connect(self):
        df = self._connect(target='tfbs', target_processors={'regulon': [to_list_nan]},
                           url='url', external_identifier='regulon', key='regulon_id',
                           explode='regulon')
        self.stack_json(df)
