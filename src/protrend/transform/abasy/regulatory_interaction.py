import pandas as pd

from protrend.io import read_from_multi_stack
from protrend.model import RegulatoryInteraction, Regulator, Gene, Organism
from protrend.transform import BaseRegulatoryInteractionTransformer
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector
from protrend.transform.abasy.gene import GeneTransformer
from protrend.transform.abasy.organism import OrganismTransformer
from protrend.transform.abasy.regulator import RegulatorTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str, regulatory_effect_abasy


class RegulatoryInteractionTransformer(AbasyTransformer, BaseRegulatoryInteractionTransformer,
                                       source='abasy',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash',
                       'id', 'regulator', 'target', 'Effect', 'Evidence', 'source', 'taxonomy',
                       'regulator_taxonomy', 'target_taxonomy'])

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network,
                                   regulator=[rstrip, lstrip],
                                   target=[rstrip, lstrip],
                                   Effect=[rstrip, lstrip])

        network = network.dropna(subset=['regulator', 'target', 'taxonomy', 'Effect'])
        network = self.drop_empty_string(network, 'regulator', 'target')
        network = self.drop_duplicates(df=network, subset=['regulator', 'target', 'taxonomy', 'Effect'],
                                       perfect_match=True)

        regulator_taxonomy = network['regulator'] + network['taxonomy']
        gene_taxonomy = network['target'] + network['taxonomy']

        network = network.assign(regulator_taxonomy=regulator_taxonomy,
                                 gene_taxonomy=gene_taxonomy)
        network = network.rename(columns={'Effect': 'regulatory_effect'})
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = self.select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism', 'ncbi_taxonomy': 'taxonomy'})
        organism = apply_processors(organism, taxonomy=to_int_str)
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = self.select_columns(regulator, 'protrend_id', 'regulator_taxonomy')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.select_columns(gene, 'protrend_id', 'gene_taxonomy')
        gene = gene.rename(columns={'protrend_id': 'gene'})
        return gene

    def transform(self) -> pd.DataFrame:
        # noinspection DuplicatedCode
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)
        organism = self.read_integrated(node='organism', columns=OrganismTransformer.columns)
        regulator = self.read_integrated(node='regulator', columns=RegulatorTransformer.columns)
        gene = self.read_integrated(node='gene', columns=GeneTransformer.columns)

        df = self._transform(network=network,
                             organism=organism, organism_key='taxonomy',
                             regulator=regulator, regulator_key='regulator_taxonomy',
                             gene=gene, gene_key='gene_taxonomy',
                             regulatory_effect_processor=regulatory_effect_abasy)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToOrganismConnector(AbasyConnector,
                                               source='abasy',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Organism,
                                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='organism')
        self.stack_json(df)


class RegulatoryInteractionToRegulatorConnector(AbasyConnector,
                                                source='abasy',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(AbasyConnector,
                                           source='abasy',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_json(df)


class RegulatorToGeneConnector(AbasyConnector,
                               source='abasy',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_json(df)
