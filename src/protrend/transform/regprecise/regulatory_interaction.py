import pandas as pd

from protrend.io import read_json_frame, read_from_stack
from protrend.model import RegulatoryInteraction, Gene, TFBS, Regulator, Effector, Organism
from protrend.transform import RegulatoryInteractionMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.effector import EffectorTransformer
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.tfbs import TFBSTransformer
from protrend.utils import SetList
from protrend.utils.processors import (apply_processors, regulatory_effect_regprecise, to_list_nan)


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, RegPreciseTransformer,
                                       source='regprecise',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=60,
                                       register=True):
    default_transform_stack = {'organism': 'integrated_organism.json',
                               'effector': 'integrated_effector.json',
                               'regulator': 'integrated_regulator.json',
                               'gene': 'integrated_gene.json',
                               'tfbs': 'integrated_tfbs.json'}
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash'])

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        network = self.select_columns(network, 'regulon_id', 'genome', 'url', 'regulation_mode', 'effector',
                                      'tfbs', 'gene')
        network = apply_processors(network, effector=to_list_nan, tfbs=to_list_nan, gene=to_list_nan)
        network = network.explode('effector')
        network = network.explode('tfbs')
        network = network.explode('gene')
        network = network.dropna(subset=['regulon_id', 'gene'])
        network = self.drop_empty_string(network, 'regulon_id', 'gene')
        network = network.rename(columns={'regulation_mode': 'regulatory_effect',
                                          'effector': 'effector_id', 'gene': 'tg_gene', 'tfbs': 'tfbs_id'})
        return network

    def transform_effector(self, effector: pd.DataFrame) -> pd.DataFrame:
        effector = self.select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'effector'})
        return effector

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = self.select_columns(regulator, 'protrend_id', 'regulon_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.select_columns(gene, 'protrend_id', 'regprecise_locus_tag')
        gene = apply_processors(gene, regprecise_locus_tag=to_list_nan)
        gene = gene.explode('regprecise_locus_tag')
        gene = gene.rename(columns={'protrend_id': 'gene', 'regprecise_locus_tag': 'tg_gene'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = self.select_columns(tfbs, 'protrend_id', 'tfbs_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform(self) -> pd.DataFrame:
        network = read_from_stack(stack=self.transform_stack, key='regulator',
                                  columns=RegulatorTransformer.read_columns, reader=read_json_frame)

        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.read_columns, reader=read_json_frame)

        effector = read_from_stack(stack=self.transform_stack, key='effector',
                                   columns=EffectorTransformer.read_columns, reader=read_json_frame)

        # noinspection DuplicatedCode
        regulator = read_from_stack(stack=self.transform_stack, key='regulator',
                                    columns=RegulatorTransformer.read_columns, reader=read_json_frame)

        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=GeneTransformer.read_columns, reader=read_json_frame)

        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs',
                               columns=TFBSTransformer.columns, reader=read_json_frame)

        df = self._transform(network=network,
                             organism=organism, organism_key='genome',
                             regulator=regulator, regulator_key='regulon_id',
                             gene=gene, gene_key='tg_gene',
                             tfbs=tfbs, tfbs_key='tfbs_id',
                             effector=effector, effector_key='effector_id',
                             regulatory_effect_processor=regulatory_effect_regprecise)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToOrganismConnector(RegPreciseConnector,
                                               source='regprecise',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Organism,
                                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='organism')
        self.stack_json(df)


class RegulatoryInteractionToEffectorConnector(RegPreciseConnector,
                                               source='regprecise',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Effector,
                                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='effector')
        self.stack_json(df)


class RegulatoryInteractionToRegulatorConnector(RegPreciseConnector,
                                                source='regprecise',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(RegPreciseConnector,
                                           source='regprecise',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_json(df)


class RegulatoryInteractionToTFBSConnector(RegPreciseConnector,
                                           source='regprecise',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='tfbs')
        self.stack_json(df)


class RegulatorToOrganismConnector(RegPreciseConnector,
                                   source='regprecise',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Organism,
                                   register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='organism')
        self.stack_json(df)


class RegulatorToEffectorConnector(RegPreciseConnector,
                                   source='regprecise',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Effector,
                                   register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='effector')
        self.stack_json(df)


class RegulatorToGeneConnector(RegPreciseConnector,
                               source='regprecise',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_json(df)


class RegulatorToTFBSConnector(RegPreciseConnector,
                               source='regprecise',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='tfbs')
        self.stack_json(df)


class GeneToOrganismConnector(RegPreciseConnector,
                              source='regprecise',
                              version='0.0.0',
                              from_node=Gene,
                              to_node=Organism,
                              register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='organism')
        self.stack_json(df)


class GeneToTFBSConnector(RegPreciseConnector,
                          source='regprecise',
                          version='0.0.0',
                          from_node=Gene,
                          to_node=TFBS,
                          register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_json(df)


class TFBSToOrganismConnector(RegPreciseConnector,
                              source='regprecise',
                              version='0.0.0',
                              from_node=TFBS,
                              to_node=Organism,
                              register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='tfbs', target_column='organism')
        self.stack_json(df)
