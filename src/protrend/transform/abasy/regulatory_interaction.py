import pandas as pd

from protrend.io import read_json_frame, read_from_stack, read_from_multi_stack
from protrend.model import RegulatoryInteraction, Regulator, Gene, Organism
from protrend.transform.abasy.base import AbasyTransformer, AbasyConnector
from protrend.transform.abasy.gene import GeneTransformer
from protrend.transform.abasy.organism import OrganismTransformer
from protrend.transform.abasy.regulator import RegulatorTransformer
from protrend.utils import SetList, build_stack
from protrend.utils.processors import apply_processors, rstrip, lstrip, to_int_str, regulatory_effect_abasy


class RegulatoryInteractionTransformer(AbasyTransformer,
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
        network = network.dropna(subset=['regulator', 'target', 'taxonomy', 'Effect'])
        network = self.drop_duplicates(df=network, subset=['regulator', 'target', 'taxonomy', 'Effect'],
                                       perfect_match=True)

        network = apply_processors(network,
                                   regulator=[rstrip, lstrip],
                                   target=[rstrip, lstrip],
                                   Effect=[rstrip, lstrip])

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
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)
        network = self.transform_network(network)

        # by organism
        organism_stack = build_stack(source=self.source, version=self.version,
                                     stack_to_load={'organism': 'integrated_organism.json'}, sa=False)
        organism = read_from_stack(stack=organism_stack,
                                   key='organism',
                                   columns=OrganismTransformer.columns,
                                   reader=read_json_frame)
        organism = self.transform_organism(organism)
        network_organism = pd.merge(network, organism, on='taxonomy')

        # by regulator
        regulator_stack = build_stack(source=self.source, version=self.version,
                                      stack_to_load={'regulator': 'integrated_regulator.json'}, sa=False)
        regulator = read_from_stack(stack=regulator_stack,
                                    key='regulator',
                                    columns=RegulatorTransformer.columns,
                                    reader=read_json_frame)
        regulator = self.transform_regulator(regulator)
        network_organism_regulator = pd.merge(network_organism, regulator, on='regulator_taxonomy')

        # by gene
        gene_stack = build_stack(source=self.source, version=self.version,
                                 stack_to_load={'gene': 'integrated_gene.json'}, sa=False)
        gene = read_from_stack(stack=gene_stack,
                               key='gene',
                               columns=GeneTransformer.columns,
                               reader=read_json_frame)
        gene = self.transform_gene(gene)
        network_organism_regulator_gene = pd.merge(network_organism_regulator, gene, on='gene_taxonomy')

        regulatory_interaction = apply_processors(network_organism_regulator_gene,
                                                  regulatory_effect=regulatory_effect_abasy)

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
