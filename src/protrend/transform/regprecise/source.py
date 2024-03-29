from typing import Dict, Callable, List

import pandas as pd

from protrend.io.utils import read_source, read_rfam
from protrend.model import (Source, Effector, Gene, Organism, Pathway, Regulator, RegulatoryFamily,
                            RegulatoryInteraction, TFBS)
from protrend.transform.mix_ins import SourceMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.regulatory_family import RegulatoryFamilyTransformer
from protrend.utils import SetList
from protrend.utils.constants import DATABASE
from protrend.utils.processors import to_list_nan, take_first


class SourceTransformer(SourceMixIn, RegPreciseTransformer,
                        source='regprecise',
                        version='0.0.0',
                        node=Source,
                        order=100,
                        register=True):
    name = ['regprecise']
    type = [DATABASE]
    url = ['https://regprecise.lbl.gov/']
    doi = ['10.1186/1471-2164-14-745']
    authors = [['Pavel S Novichkov', 'Alexey E Kazakov', 'Dmitry A Ravcheev', 'Semen A Leyn', 'Galina Y Kovaleva',
                'Roman A Sutormin', 'Marat D Kazanov', 'William Riehl', 'Adam P Arkin',
                'Inna Dubchak', 'Dmitry A Rodionov']]
    description = ['RegPrecise 3.0: A resource for genome-scale exploration of transcriptional regulation in bacteria']

    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])


class SourceConnector(RegPreciseConnector, register=False):

    def _connect(self,
                 target: str,
                 target_processors: Dict[str, List[Callable]],
                 url: str,
                 external_identifier: str,
                 key: str):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target=target,
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors=target_processors)

        source_ids, target_ids = self.merge_source_target(source_df=source_df, target_df=target_df,
                                                          cardinality='one_to_many')

        if url in target_df.columns and external_identifier in target_df.columns:
            size = len(target_ids)
            kwargs = dict(url=target_df[url].to_list(),
                          external_identifier=target_df[external_identifier].to_list(),
                          key=[key] * size)

        else:
            kwargs = {}

        return self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)


class SourceToEffectorConnector(SourceConnector,
                                source='regprecise',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Effector,
                                register=True):

    def connect(self):
        df = self._connect(target='effector', target_processors={},
                           url='url', external_identifier='effector_id', key='effector_id')
        self.stack_connections(df)


class SourceToGeneConnector(SourceConnector,
                            source='regprecise',
                            version='0.0.0',
                            from_node=Source,
                            to_node=Gene,
                            register=True):

    def connect(self):
        source_df, target_df = self.transform_stacks(source='source',
                                                     target='gene',
                                                     source_column='protrend_id',
                                                     target_column='protrend_id',
                                                     source_processors={},
                                                     target_processors={'regulon': [to_list_nan]})

        if 'regulon' in target_df:
            target_df = target_df.explode('regulon')

        source_df = source_df.reset_index(drop=True)
        target_df = target_df.reset_index(drop=True)
        df = pd.concat([source_df, target_df], axis=1)

        target_ids = df['target_col'].to_list()

        source_ids = df['source_col'].dropna().to_list()
        source_ids *= len(target_ids)

        if 'url' in target_df.columns and 'regulon' in target_df.columns:
            size = len(target_ids)
            kwargs = dict(url=target_df['url'].to_list(),
                          external_identifier=target_df['regulon'].to_list(),
                          key=['regulon_id'] * size)

        else:
            kwargs = {}

        df = self.connection_frame(source_ids=source_ids, target_ids=target_ids, kwargs=kwargs)
        self.stack_connections(df)


class SourceToOrganismConnector(SourceConnector,
                                source='regprecise',
                                version='0.0.0',
                                from_node=Source,
                                to_node=Organism,
                                register=True):

    def connect(self):
        df = self._connect(target='organism', target_processors={},
                           url='url', external_identifier='genome_id', key='genome_id')
        self.stack_connections(df)


class SourceToPathwayConnector(SourceConnector,
                               source='regprecise',
                               version='0.0.0',
                               from_node=Source,
                               to_node=Pathway,
                               register=True):

    def connect(self):
        df = self._connect(target='pathway', target_processors={},
                           url='url', external_identifier='pathway_id', key='pathway_id')
        self.stack_connections(df)


class SourceToRegulatorConnector(SourceConnector,
                                 source='regprecise',
                                 version='0.0.0',
                                 from_node=Source,
                                 to_node=Regulator,
                                 register=True):

    def connect(self):
        df = self._connect(target='regulator', target_processors={},
                           url='url', external_identifier='regulon_id', key='regulon_id')
        self.stack_connections(df)


class SourceToRegulatoryFamilyConnector(RegPreciseConnector,
                                        source='regprecise',
                                        version='0.0.0',
                                        from_node=Source,
                                        to_node=RegulatoryFamily,
                                        register=True):

    def connect(self):
        source = read_source(source=self.source, version=self.version, columns=SourceTransformer.columns)
        target = read_rfam(source=self.source, version=self.version, columns=RegulatoryFamilyTransformer.columns)

        protrend_id = source.loc[0, 'protrend_id']

        from_identifiers = []
        to_identifiers = []
        urls = []
        external_ids = []
        keys = []
        for target_key in ('tffamily_id', 'collection_id', 'riboswitch_id'):
            target_df = target.dropna(subset=[target_key])

            to_identifier = target_df['protrend_id'].to_list()

            size = len(to_identifier)
            from_identifier = [protrend_id] * size


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
        self.stack_connections(df)


class SourceToRegulatoryInteractionConnector(SourceConnector,
                                             source='regprecise',
                                             version='0.0.0',
                                             from_node=Source,
                                             to_node=RegulatoryInteraction,
                                             register=True):

    def connect(self):
        df = self._connect(target='rin', target_processors={},
                           url='url', external_identifier='regulon_id', key='regulon_id')
        self.stack_connections(df)


class SourceToTFBSConnector(SourceConnector,
                            source='regprecise',
                            version='0.0.0',
                            from_node=Source,
                            to_node=TFBS,
                            register=True):

    def connect(self):
        df = self._connect(target='tfbs', target_processors={'regulon': [to_list_nan, take_first]},
                           url='url', external_identifier='regulon', key='regulon_id')
        self.stack_connections(df)
