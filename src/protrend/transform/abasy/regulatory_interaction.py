import pandas as pd

from protrend.io.utils import read_organism, read_regulator, read_gene
from protrend.model import RegulatoryInteraction, Regulator, Gene, Organism
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector, read_abasy_networks
from protrend.transform.abasy.gene import GeneTransformer
from protrend.transform.abasy.organism import OrganismTransformer
from protrend.transform.abasy.regulator import RegulatorTransformer
from protrend.transform.mix_ins import RegulatoryInteractionMixIn
from protrend.transform.transformations import select_columns, drop_empty_string, drop_duplicates
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str, regulatory_effect_abasy


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, AbasyTransformer,
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
        network = drop_empty_string(network, 'regulator', 'target')
        network = drop_duplicates(df=network, subset=['regulator', 'target', 'taxonomy', 'Effect'],
                                  perfect_match=True)

        regulator_taxonomy = network['regulator'] + network['taxonomy']
        gene_taxonomy = network['target'] + network['taxonomy']

        network = network.assign(regulator_taxonomy=regulator_taxonomy,
                                 gene_taxonomy=gene_taxonomy)
        network = network.rename(columns={'Effect': 'regulatory_effect'})
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism', 'ncbi_taxonomy': 'taxonomy'})
        organism = apply_processors(organism, taxonomy=to_int_str)
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'protrend_id', 'regulator_taxonomy')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'gene_taxonomy')
        gene = gene.rename(columns={'protrend_id': 'gene'})
        return gene

    def transform(self) -> pd.DataFrame:
        network = read_abasy_networks(self.source, self.version)

        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)

        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)

        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)

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
        self.stack_connections(df)


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
        self.stack_connections(df)


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
        self.stack_connections(df)


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
        self.stack_connections(df)
