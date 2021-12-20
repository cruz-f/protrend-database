import pandas as pd

from protrend.io import read_json_frame, read_json_lines, read_from_stack
from protrend.model import RegulatoryInteraction, Regulator, Gene, TFBS
from protrend.transform import BaseRegulatoryInteractionTransformer
from protrend.transform.collectf.base import CollectfTransformer, CollectfConnector
from protrend.transform.collectf.gene import GeneTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.transform.collectf.tfbs import TFBSTransformer
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, regulatory_effect_collectf, to_list_nan


class RegulatoryInteractionTransformer(CollectfTransformer, BaseRegulatoryInteractionTransformer,
                                       source='collectf',
                                       version='0.0.1',
                                       node=RegulatoryInteraction,
                                       order=70,
                                       register=True):
    default_transform_stack = {'regulator': 'integrated_regulator.json',
                               'gene': 'integrated_gene.json',
                               'tfbs': 'integrated_tfbs.json',
                               'network': 'TFBS.json'}
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash'])

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        network = self.select_columns(network, 'tfbs_id', 'regulon', 'gene', 'mode')
        network = apply_processors(network, regulon=to_list_nan, gene=to_list_nan)
        network = network.explode('regulon')
        network = network.explode('gene')
        network = network.dropna(subset=['regulon', 'gene'])
        network = network.rename(columns={'mode': 'regulatory_effect', 'gene': 'tg_gene'})
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = self.select_columns(organism, 'organism_protrend_id', 'uniprot_accession')
        organism = organism.rename(columns={'organism_protrend_id': 'organism', 'uniprot_accession': 'regulon'})
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = self.select_columns(regulator, 'protrend_id', 'uniprot_accession')
        regulator = regulator.rename(columns={'protrend_id': 'regulator', 'uniprot_accession': 'regulon'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.select_columns(gene, 'protrend_id', 'locus_tag_old')
        gene = gene.rename(columns={'protrend_id': 'gene', 'locus_tag_old': 'tg_gene'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = self.select_columns(tfbs, 'tfbs_id', 'protrend_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform(self) -> pd.DataFrame:
        network = read_from_stack(stack=self.transform_stack, key='network',
                                  columns=TFBSTransformer.read_columns, reader=read_json_lines)

        organism = read_from_stack(stack=self.transform_stack, key='regulator',
                                   columns=RegulatorTransformer.read_columns, reader=read_json_frame)

        regulator = read_from_stack(stack=self.transform_stack, key='regulator',
                                    columns=RegulatorTransformer.read_columns, reader=read_json_frame)

        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=GeneTransformer.read_columns, reader=read_json_frame)

        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs',
                               columns=TFBSTransformer.columns, reader=read_json_frame)

        df = self._transform(network=network,
                             organism=organism, organism_key='regulon',
                             regulator=regulator, regulator_key='regulon',
                             gene=gene, gene_key='tg_gene',
                             tfbs=tfbs, tfbs_key='tfbs_id',
                             regulatory_effect_processor=regulatory_effect_collectf)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToRegulatorConnector(CollectfConnector,
                                                source='collectf',
                                                version='0.0.1',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(CollectfConnector,
                                           source='collectf',
                                           version='0.0.1',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_json(df)


class RegulatoryInteractionToTFBSConnector(CollectfConnector,
                                           source='collectf',
                                           version='0.0.1',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='tfbs')
        self.stack_json(df)


class RegulatorToGeneConnector(CollectfConnector,
                               source='collectf',
                               version='0.0.1',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_json(df)


class RegulatorToTFBSConnector(CollectfConnector,
                               source='collectf',
                               version='0.0.1',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='tfbs')
        self.stack_json(df)


class GeneToTFBSConnector(CollectfConnector,
                          source='collectf',
                          version='0.0.1',
                          from_node=Gene,
                          to_node=TFBS,
                          register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_json(df)
