import pandas as pd

from protrend.io.utils import read_regulator, read_organism, read_gene, read_tfbs, read_effector
from protrend.model import RegulatoryInteraction, Gene, TFBS, Regulator, Effector, Organism
from protrend.transform.mix_ins import RegulatoryInteractionMixIn
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.effector import EffectorTransformer
from protrend.transform.regprecise.gene import GeneTransformer
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.tfbs import TFBSTransformer
from protrend.transform.transformations import select_columns, drop_empty_string
from protrend.utils import SetList
from protrend.utils.processors import (apply_processors, regulatory_effect_regprecise, to_list_nan)


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, RegPreciseTransformer,
                                       source='regprecise',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=60,
                                       register=True):
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash'])

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        network = select_columns(network, 'regulon_id', 'genome', 'url', 'regulation_mode', 'effector',
                                 'tfbs', 'gene')
        network = apply_processors(network, effector=to_list_nan, tfbs=to_list_nan, gene=to_list_nan)
        network = network.explode('effector')
        network = network.explode('tfbs')
        network = network.explode('gene')
        network = network.dropna(subset=['regulon_id', 'gene'])
        network = drop_empty_string(network, 'regulon_id', 'gene')
        network = network.rename(columns={'regulation_mode': 'regulatory_effect',
                                          'effector': 'effector_id', 'gene': 'tg_gene', 'tfbs': 'tfbs_id'})
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id', 'genome_id')
        organism = organism.rename(columns={'protrend_id': 'organism', 'genome_id': 'genome'})
        return organism

    def transform_effector(self, effector: pd.DataFrame) -> pd.DataFrame:
        effector = select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'effector'})
        return effector

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'protrend_id', 'regulon_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'regprecise_locus_tag')
        gene = apply_processors(gene, regprecise_locus_tag=to_list_nan)
        gene = gene.explode('regprecise_locus_tag')
        gene = gene.rename(columns={'protrend_id': 'gene', 'regprecise_locus_tag': 'tg_gene'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = select_columns(tfbs, 'protrend_id', 'tfbs_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform(self) -> pd.DataFrame:
        network = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)
        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)
        regulator = network.copy()
        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)
        tfbs = read_tfbs(source=self.source, version=self.version, columns=TFBSTransformer.columns)
        effector = read_effector(source=self.source, version=self.version, columns=EffectorTransformer.columns)

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

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='organism')
        self.stack_connections(df)


class RegulatoryInteractionToEffectorConnector(RegPreciseConnector,
                                               source='regprecise',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Effector,
                                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='effector')
        self.stack_connections(df)


class RegulatoryInteractionToRegulatorConnector(RegPreciseConnector,
                                                source='regprecise',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_connections(df)


class RegulatoryInteractionToGeneConnector(RegPreciseConnector,
                                           source='regprecise',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_connections(df)


class RegulatoryInteractionToTFBSConnector(RegPreciseConnector,
                                           source='regprecise',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='tfbs')
        self.stack_connections(df)


class RegulatorToOrganismConnector(RegPreciseConnector,
                                   source='regprecise',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Organism,
                                   register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='organism')
        self.stack_connections(df)


class RegulatorToEffectorConnector(RegPreciseConnector,
                                   source='regprecise',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Effector,
                                   register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='effector')
        self.stack_connections(df)


class RegulatorToGeneConnector(RegPreciseConnector,
                               source='regprecise',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_connections(df)


class RegulatorToTFBSConnector(RegPreciseConnector,
                               source='regprecise',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='tfbs')
        self.stack_connections(df)


class GeneToOrganismConnector(RegPreciseConnector,
                              source='regprecise',
                              version='0.0.0',
                              from_node=Gene,
                              to_node=Organism,
                              register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='organism')
        self.stack_connections(df)


class GeneToTFBSConnector(RegPreciseConnector,
                          source='regprecise',
                          version='0.0.0',
                          from_node=Gene,
                          to_node=TFBS,
                          register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_connections(df)


class TFBSToOrganismConnector(RegPreciseConnector,
                              source='regprecise',
                              version='0.0.0',
                              from_node=TFBS,
                              to_node=Organism,
                              register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='tfbs', target_column='organism')
        self.stack_connections(df)
