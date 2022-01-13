import pandas as pd

from protrend.io import read_json_lines, read
from protrend.io.utils import read_regulator, read_gene, read_tfbs
from protrend.model import RegulatoryInteraction, Regulator, Gene, TFBS
from protrend.transform.collectf.base import CollecTFTransformer, CollecTFConnector
from protrend.transform.collectf.gene import GeneTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.transform.mix_ins import RegulatoryInteractionMixIn
from protrend.transform.transformations import select_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, regulatory_effect_collectf, to_list_nan


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, CollecTFTransformer,
                                       source='collectf',
                                       version='0.0.1',
                                       node=RegulatoryInteraction,
                                       order=70,
                                       register=True):
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash'])

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        network = select_columns(network, 'tfbs_id', 'regulon', 'gene', 'mode')
        network = apply_processors(network, regulon=to_list_nan, gene=to_list_nan)
        network = network.explode('regulon')
        network = network.explode('gene')
        network = network.dropna(subset=['regulon', 'gene'])
        network = network.rename(columns={'mode': 'regulatory_effect', 'gene': 'tg_gene'})
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'organism_protrend_id', 'uniprot_accession')
        organism = organism.rename(columns={'organism_protrend_id': 'organism', 'uniprot_accession': 'regulon'})
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'protrend_id', 'uniprot_accession')
        regulator = regulator.rename(columns={'protrend_id': 'regulator', 'uniprot_accession': 'regulon'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'locus_tag_old')
        gene = gene.rename(columns={'protrend_id': 'gene', 'locus_tag_old': 'tg_gene'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = select_columns(tfbs, 'tfbs_id', 'protrend_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform(self) -> pd.DataFrame:
        network = read(source=self.source, version=self.version, file='TFBS.json',
                       reader=read_json_lines, default=pd.DataFrame(columns=TFBSTransformer.columns))

        organism = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)

        regulator = organism.copy()

        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)

        tfbs = read_tfbs(source=self.source, version=self.version, columns=TFBSTransformer.columns)

        df = self._transform(network=network,
                             organism=organism, organism_key='regulon',
                             regulator=regulator, regulator_key='regulon',
                             gene=gene, gene_key='tg_gene',
                             tfbs=tfbs, tfbs_key='tfbs_id',
                             regulatory_effect_processor=regulatory_effect_collectf)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToRegulatorConnector(CollecTFConnector,
                                                source='collectf',
                                                version='0.0.1',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_connections(df)


class RegulatoryInteractionToGeneConnector(CollecTFConnector,
                                           source='collectf',
                                           version='0.0.1',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_connections(df)


class RegulatoryInteractionToTFBSConnector(CollecTFConnector,
                                           source='collectf',
                                           version='0.0.1',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='tfbs')
        self.stack_connections(df)


class RegulatorToGeneConnector(CollecTFConnector,
                               source='collectf',
                               version='0.0.1',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_connections(df)


class RegulatorToTFBSConnector(CollecTFConnector,
                               source='collectf',
                               version='0.0.1',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='tfbs')
        self.stack_connections(df)


class GeneToTFBSConnector(CollecTFConnector,
                          source='collectf',
                          version='0.0.1',
                          from_node=Gene,
                          to_node=TFBS,
                          register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_connections(df)
