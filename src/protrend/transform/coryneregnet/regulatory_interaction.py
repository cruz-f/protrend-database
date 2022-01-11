import pandas as pd

from protrend.io import read_from_multi_stack
from protrend.model import RegulatoryInteraction, Regulator, Gene, TFBS
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer, CoryneRegNetConnector
from protrend.transform.coryneregnet.gene import GeneTransformer
from protrend.transform.coryneregnet.organism import OrganismTransformer
from protrend.transform.coryneregnet.regulator import RegulatorTransformer
from protrend.transform.coryneregnet.tfbs import TFBSTransformer
from protrend.transform.mix_ins import RegulatoryInteractionMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, select_columns
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, regulatory_effect_coryneregnet, rstrip, lstrip, to_int_str


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, CoryneRegNetTransformer,
                                       source='coryneregnet',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence',
                       'PMID', 'Source', 'taxonomy', 'source'])

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        network = apply_processors(network,
                                   TF_locusTag=[rstrip, lstrip],
                                   TG_locusTag=[rstrip, lstrip],
                                   Role=[rstrip, lstrip])

        network = network.dropna(subset=['TF_locusTag', 'TG_locusTag'])
        network = drop_empty_string(network, 'TF_locusTag', 'TG_locusTag')
        network = drop_duplicates(df=network,
                                       subset=['TF_locusTag', 'TG_locusTag', 'Binding_site', 'Role'],
                                       perfect_match=True)

        tfbs_id = network['TF_locusTag'] + network['TG_locusTag'] + network['Binding_site'] + network['taxonomy']
        network = network.assign(tfbs_id=tfbs_id)

        network = network.rename(columns={'Role': 'regulatory_effect'})
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = organism.rename(columns={'protrend_id': 'organism', 'ncbi_taxonomy': 'taxonomy'})
        organism = apply_processors(organism, taxonomy=to_int_str)
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'protrend_id', 'TF_locusTag')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'TG_locusTag')
        gene = gene.rename(columns={'protrend_id': 'gene'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = select_columns(tfbs, 'protrend_id', 'TF_locusTag', 'TG_locusTag', 'Binding_site', 'taxonomy')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})

        tfbs_id = tfbs['TF_locusTag'] + tfbs['TG_locusTag'] + tfbs['Binding_site'] + tfbs['taxonomy']
        tfbs = tfbs.assign(tfbs_id=tfbs_id)
        return tfbs

    def transform(self) -> pd.DataFrame:
        # noinspection DuplicatedCode
        network = read_from_multi_stack(stack=self.transform_stack, key='network', columns=self.default_network_columns)
        organism = self.read_integrated(node='organism', columns=OrganismTransformer.columns)
        regulator = self.read_integrated(node='regulator', columns=RegulatorTransformer.columns)
        gene = self.read_integrated(node='gene', columns=GeneTransformer.columns)
        tfbs = self.read_integrated(node='tfbs', columns=TFBSTransformer.columns)

        df = self._transform(network=network,
                             organism=organism, organism_key='taxonomy',
                             regulator=regulator, regulator_key='TF_locusTag',
                             gene=gene, gene_key='TG_locusTag',
                             tfbs=tfbs, tfbs_key='tfbs_id',
                             regulatory_effect_processor=regulatory_effect_coryneregnet)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToRegulatorConnector(CoryneRegNetConnector,
                                                source='coryneregnet',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(CoryneRegNetConnector,
                                           source='coryneregnet',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_json(df)


class RegulatoryInteractionToTFBSConnector(CoryneRegNetConnector,
                                           source='coryneregnet',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='tfbs')
        self.stack_json(df)


class RegulatorToGeneConnector(CoryneRegNetConnector,
                               source='coryneregnet',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_json(df)


class RegulatorToTFBSConnector(CoryneRegNetConnector,
                               source='coryneregnet',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='tfbs')
        self.stack_json(df)


class GeneToTFBSConnector(CoryneRegNetConnector,
                          source='coryneregnet',
                          version='0.0.0',
                          from_node=Gene,
                          to_node=TFBS,
                          register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_json(df)
