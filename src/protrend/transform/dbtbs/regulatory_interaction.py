import pandas as pd

from protrend.io import read_json_lines, read
from protrend.io.utils import read_organism, read_regulator, read_gene, read_tfbs
from protrend.model import RegulatoryInteraction, Regulator, Gene, TFBS
from protrend.transform.dbtbs.base import DBTBSTransformer, DBTBSConnector
from protrend.transform.dbtbs.gene import GeneTransformer
from protrend.transform.dbtbs.organism import OrganismTransformer
from protrend.transform.dbtbs.regulator import RegulatorTransformer
from protrend.transform.dbtbs.tfbs import TFBSTransformer
from protrend.transform.mix_ins import RegulatoryInteractionMixIn
from protrend.transform.transformations import select_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, regulatory_effect_dbtbs, to_list_nan


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, DBTBSTransformer,
                                       source='dbtbs',
                                       version='0.0.4',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):

    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash'])

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        network = select_columns(network, 'identifier', 'url', 'regulation', 'tf', 'gene')
        network = apply_processors(network,
                                   url=to_list_nan, regulation=to_list_nan, tf=to_list_nan, gene=to_list_nan)
        network = network.explode('url')
        network = network.explode('regulation')
        network = network.explode('tf')
        network = network.explode('gene')
        network = network.dropna(subset=['tf', 'gene'])
        network = network.rename(columns={'regulation': 'regulatory_effect', 'gene': 'tg_gene'})
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id', 'identifier')
        organism = organism.rename(columns={'protrend_id': 'organism'})
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'protrend_id', 'dbtbs_name')
        regulator = regulator.rename(columns={'protrend_id': 'regulator', 'dbtbs_name': 'tf'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'dbtbs_name')
        gene = gene.rename(columns={'protrend_id': 'gene', 'dbtbs_name': 'tg_gene'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = select_columns(tfbs, 'identifier', 'protrend_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform(self) -> pd.DataFrame:
        network = read(source=self.source, version=self.version, file='TFBS.json',
                       reader=read_json_lines, default=pd.DataFrame(columns=TFBSTransformer.columns))
        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)
        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)
        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)
        tfbs = read_tfbs(source=self.source, version=self.version, columns=TFBSTransformer.columns)

        network_ids = network['identifier'].to_list()

        organism_id = organism.loc[0, 'protrend_id']
        organism_ids = [organism_id] * len(network_ids)

        organism = pd.DataFrame({'protrend_id': organism_ids, 'identifier': network_ids})

        df = self._transform(network=network,
                             organism=organism, organism_key='identifier',
                             regulator=regulator, regulator_key='tf',
                             gene=gene, gene_key='tg_gene',
                             tfbs=tfbs, tfbs_key='identifier',
                             regulatory_effect_processor=regulatory_effect_dbtbs)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToRegulatorConnector(DBTBSConnector,
                                                source='dbtbs',
                                                version='0.0.4',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_connections(df)


class RegulatoryInteractionToGeneConnector(DBTBSConnector,
                                           source='dbtbs',
                                           version='0.0.4',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_connections(df)


class RegulatoryInteractionToTFBSConnector(DBTBSConnector,
                                           source='dbtbs',
                                           version='0.0.4',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='tfbs')
        self.stack_connections(df)


class RegulatorToGeneConnector(DBTBSConnector,
                               source='dbtbs',
                               version='0.0.4',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_connections(df)


class RegulatorToTFBSConnector(DBTBSConnector,
                               source='dbtbs',
                               version='0.0.4',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='tfbs')
        self.stack_connections(df)


class GeneToTFBSConnector(DBTBSConnector,
                          source='dbtbs',
                          version='0.0.4',
                          from_node=Gene,
                          to_node=TFBS,
                          register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_connections(df)
