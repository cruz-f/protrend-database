import pandas as pd

from protrend import Organism
from protrend.io import read_json_frame, read_from_stack
from protrend.model import RegulatoryInteraction, Regulator, Gene
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector, read_abasy_network
from protrend.transform.abasy.gene import GeneTransformer
from protrend.transform.abasy.organism import OrganismTransformer
from protrend.transform.abasy.regulator import RegulatorTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str, regulatory_effect_abasy


class RegulatoryInteractionTransformer(AbasyTransformer,
                                       source='abasy',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'gene': 'integrated_gene.json',
                               'organism': 'integrated_organism.json'}
    columns = SetList(['regulator_effector', 'regulator', 'operon', 'genes', 'tfbss', 'regulatory_effect',
                       'regulatory_interaction_hash', 'protrend_id',
                       'id', 'source', 'target', 'Effect', 'Evidence', 'taxonomy', 'source_target_taxonomy'])

    def transform_networks(self, networks: pd.DataFrame) -> pd.DataFrame:
        networks = networks.dropna(subset=['source', 'target', 'taxonomy', 'Effect'])
        networks = self.drop_duplicates(df=networks, subset=['source', 'target', 'taxonomy', 'Effect'],
                                        perfect_match=True)

        networks = apply_processors(networks,
                                    source=[rstrip, lstrip],
                                    target=[rstrip, lstrip],
                                    Effect=[rstrip, lstrip])

        regulator_taxonomy = networks['source'] + networks['taxonomy']
        gene_taxonomy = networks['target'] + networks['taxonomy']

        networks = networks.assign(regulator_taxonomy=regulator_taxonomy,
                                   gene_taxonomy=gene_taxonomy)

        return networks

    def transform_organism(self) -> pd.DataFrame:
        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.columns, reader=read_json_frame)
        organism = self.select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism', 'ncbi_taxonomy': 'taxonomy'})
        organism = apply_processors(organism, taxonomy=to_int_str)
        return organism

    def transform_regulator(self) -> pd.DataFrame:
        regulator = read_from_stack(stack=self.transform_stack, key='regulator',
                                    columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'regulator_taxonomy')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self) -> pd.DataFrame:
        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'protrend_id', 'gene_taxonomy')
        gene = gene.rename(columns={'protrend_id': 'gene'})
        return gene

    def transform(self) -> pd.DataFrame:
        networks = self.contact_stacks(stack=self.network_stack,
                                       taxa=self.taxa_to_organism_code,
                                       default_columns=self.default_network_columns,
                                       reader=read_abasy_network)
        network = self.transform_networks(networks)

        organism = self.transform_organism()
        network_organism = pd.merge(network, organism, on='taxonomy')

        regulator = self.transform_regulator()
        network_organism_regulator = pd.merge(network_organism, regulator, on='regulator_taxonomy')

        gene = self.transform_gene()
        network_organism_regulator_gene = pd.merge(network_organism_regulator, gene, on='gene_taxonomy')

        regulatory_interaction = network_organism_regulator_gene.rename(columns={'Effect': 'regulatory_effect'})
        regulatory_interaction = apply_processors(regulatory_interaction, regulatory_effect=regulatory_effect_abasy)

        regulatory_interaction = regulatory_interaction.assign(tfbs=None, effector=None)

        regulatory_interaction = self.interaction_hash(regulatory_interaction)

        self.stack_transformed_nodes(regulatory_interaction)

        return regulatory_interaction


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
