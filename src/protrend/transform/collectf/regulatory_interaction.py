import pandas as pd

from protrend.io import read_json_lines, read
from protrend.io.utils import read_regulator, read_gene, read_tfbs, read_organism
from protrend.model import RegulatoryInteraction, Regulator, Gene, TFBS
from protrend.transform.collectf.base import CollecTFTransformer, CollecTFConnector
from protrend.transform.collectf.organism import OrganismTransformer
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
        # tfbs_id, site_start, site_end, site_strand, mode, sequence, pubmed, organism, regulon, operon, gene,
        # experimental_evidence

        network = select_columns(network, 'tfbs_id', 'mode', 'pubmed', 'organism', 'regulon',
                                 'gene', 'experimental_evidence')
        network = apply_processors(network,
                                   regulon=to_list_nan,
                                   gene=to_list_nan,
                                   pubmed=to_list_nan,
                                   experimental_evidence=to_list_nan)
        network = network.explode('regulon')
        network = network.explode('gene')
        network = network.dropna(subset=['regulon', 'gene'])
        network = network.rename(columns={'mode': 'regulatory_effect',
                                          'organism': 'organism_name',
                                          'gene': 'gene_id'})
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id', 'name')
        organism = organism.rename(columns={'protrend_id': 'organism', 'name': 'organism_name'})
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'protrend_id', 'tfbs')
        regulator = apply_processors(regulator,
                                     tfbs=to_list_nan)
        regulator = regulator.explode('tfbs')
        regulator = regulator.dropna(subset=['protrend_id', 'tfbs'])
        regulator = regulator.rename(columns={'protrend_id': 'regulator',
                                              'tfbs': 'tfbs_id'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'gene')
        gene = gene.rename(columns={'protrend_id': 'gene',
                                    'gene': 'gene_id'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = select_columns(tfbs, 'protrend_id', 'tfbs_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform(self) -> pd.DataFrame:
        # tfbs_id, site_start, site_end, site_strand, mode, sequence, pubmed, organism, regulon, operon, gene,
        # experimental_evidence
        network = read(source=self.source, version=self.version, file='TFBS.json',
                       reader=read_json_lines, default=pd.DataFrame(columns=TFBSTransformer.columns))

        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)

        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)

        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)

        tfbs = read_tfbs(source=self.source, version=self.version, columns=TFBSTransformer.columns)

        df = self._transform(network=network,
                             organism=organism, organism_key='organism_name',
                             regulator=regulator, regulator_key='tfbs_id',
                             gene=gene, gene_key='gene_id',
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

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_connections(df)
