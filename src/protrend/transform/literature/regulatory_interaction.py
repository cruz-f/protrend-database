import pandas as pd

from protrend.io.utils import read_organism, read_regulator, read_gene, read_effector
from protrend.model import RegulatoryInteraction, Effector, Regulator, Gene
from protrend.transform.literature.base import LiteratureTransformer, LiteratureConnector, read_literature_networks
from protrend.transform.literature.effector import EffectorTransformer
from protrend.transform.literature.gene import GeneTransformer
from protrend.transform.literature.organism import OrganismTransformer
from protrend.transform.literature.regulator import RegulatorTransformer
from protrend.transform.mix_ins import RegulatoryInteractionMixIn
from protrend.transform.transformations import select_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, regulatory_effect_literature, to_int_str


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, LiteratureTransformer,
                                       source='literature',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=90,
                                       register=True):
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash',
                       'regulator_locus_tag', 'gene_locus_tag',
                       'regulatory_effect', 'effector', 'mechanism',
                       'taxonomy', 'source'])

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism', 'ncbi_taxonomy': 'taxonomy'})
        organism = apply_processors(organism, taxonomy=to_int_str)
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'protrend_id', 'regulator_locus_tag')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'gene_locus_tag')
        gene = gene.rename(columns={'protrend_id': 'gene'})
        return gene

    def transform_effector(self, effector: pd.DataFrame) -> pd.DataFrame:
        effector = select_columns(effector, 'protrend_id', 'effector_name')
        effector = effector.rename(columns={'protrend_id': 'effector'})
        return effector

    def transform(self) -> pd.DataFrame:
        network = read_literature_networks(source=self.source, version=self.version)
        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)
        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)
        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)
        effector = read_effector(source=self.source, version=self.version, columns=EffectorTransformer.columns)

        df = self._transform(network=network,
                             organism=organism, organism_key='taxonomy',
                             regulator=regulator, regulator_key='regulator_locus_tag',
                             gene=gene, gene_key='gene_locus_tag',
                             effector=effector, effector_key='effector_name',
                             regulatory_effect_processor=regulatory_effect_literature)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToEffectorConnector(LiteratureConnector,
                                               source='literature',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Effector,
                                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='effector')
        self.stack_connections(df)


class RegulatoryInteractionToRegulatorConnector(LiteratureConnector,
                                                source='literature',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_connections(df)


class RegulatoryInteractionToGeneConnector(LiteratureConnector,
                                           source='literature',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_connections(df)


class RegulatorToEffectorConnector(LiteratureConnector,
                                   source='literature',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Effector,
                                   register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='effector')
        self.stack_connections(df)


class RegulatorToGeneConnector(LiteratureConnector,
                               source='literature',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_connections(df)
