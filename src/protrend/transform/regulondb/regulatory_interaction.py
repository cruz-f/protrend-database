from typing import List

import pandas as pd

from protrend.io import read
from protrend.io.utils import read_organism, read_regulator, read_gene, read_tfbs, read_effector
from protrend.model import RegulatoryInteraction, Regulator, TFBS, Gene, Effector
from protrend.transform.mix_ins import RegulatoryInteractionMixIn
from protrend.transform.regulondb.base import RegulonDBTransformer, RegulonDBConnector, regulondb_reader
from protrend.transform.regulondb.effector import EffectorTransformer
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.regulondb.organism import OrganismTransformer
from protrend.transform.regulondb.regulator import RegulatorTransformer
from protrend.transform.regulondb.tfbs import TFBSTransformer
from protrend.transform.transformations import select_columns, drop_empty_string, drop_duplicates
from protrend.utils import SetList
from protrend.utils.processors import (apply_processors, regulatory_effect_regulondb, to_int_str)


def transform_regulon_db_dataset(df: pd.DataFrame,
                                 selection: List[str],
                                 duplicates: List[str],
                                 nan: List[str]) -> pd.DataFrame:
    df = select_columns(df, *selection)
    df = df.dropna(subset=nan)
    df = drop_empty_string(df, *nan)
    df = drop_duplicates(df, subset=duplicates, perfect_match=True)
    return df


REGULONDB_STACK = dict(genetic_network=('genetic_network.txt', 43,
                                        ['regulator_id', 'regulator_name', 'regulated_id', 'regulated_name',
                                         'function_interaction', 'evidence', 'regulator_type', 'regulated_type']),
                       regulatory_interaction=('regulatory_interaction.txt', 42,
                                               ['regulatory_interaction_id', 'conformation_id', 'promoter_id',
                                                'site_id', 'ri_function', 'center_position', 'ri_dist_first_gene',
                                                'ri_first_gene_id', 'affinity_exp', 'regulatory_interaction_note',
                                                'ri_internal_comment', 'key_id_org', 'ri_sequence', 'ri_orientation',
                                                'ri_sequence_orientation']),
                       srna_interaction=('srna_interaction.txt', 38,
                                         ['srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id',
                                          'srna_function', 'srna_posleft', 'srna_posright', 'srna_sequence',
                                          'srna_regulation_type', 'srna_mechanis', 'srna_note']),
                       tf_gene=('tf_gene_interaction.txt', 43,
                                ['regulatory_interaction_id', 'conformation_id', 'object_id', 'site_id',
                                 'ri_function', 'center_position', 'ri_dist_first_gene', 'ri_first_gene_id',
                                 'affinity_exp', 'regulatory_interaction_note', 'ri_internal_comment', 'key_id_org',
                                 'ri_sequence', 'ri_orientation', 'ri_sequence_orientation', 'object_type']),
                       tf=('transcription_factor.txt', 38,
                           ['transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                            'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                            'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence']),
                       sigma=('sigma_tmp.txt', 36,
                              ['sigma_id', 'sigma_name', 'sigma_synonyms', 'sigma_gene_id', 'sigma_gene_name',
                               'sigma_coregulators', 'sigma_notes', 'sigma_sigmulon_genes', 'key_id_org']),
                       srna=('srna_interaction.txt', 38,
                             ['srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id',
                              'srna_function', 'srna_posleft', 'srna_posright', 'srna_sequence',
                              'srna_regulation_type', 'srna_mechanis', 'srna_note']),
                       target_gene=('gene.txt', 39,
                                    ['gene_id', 'gene_name', 'gene_posleft', 'gene_posright',
                                     'gene_strand', 'gene_sequence', 'gc_content', 'cri_score',
                                     'gene_note', 'gene_internal_comment', 'key_id_org', 'gene_type']),
                       site=('site.txt', 35, ['site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
                                              'site_internal_comment', 'key_id_org', 'site_length']),
                       regulator_effector=('effector.txt', 34,
                                           ['effector_id', 'effector_name', 'category', 'effector_type',
                                            'effector_note', 'effector_internal_comment', 'key_id_org']),
                       conformation=('conformation.txt', 36,
                                     ['conformation_id', 'transcription_factor_id', 'final_state', 'conformation_note',
                                      'interaction_type', 'conformation_internal_comment', 'key_id_org',
                                      'conformation_type', 'apo_holo_conformation']),
                       conformation_effector=('conformation_effector_link.txt', 29,
                                              ['effector_id', 'conformation_id']))


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, RegulonDBTransformer,
                                       source='regulondb',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):
    columns = SetList(['protrend_id', 'organism', 'regulator', 'gene', 'tfbs', 'effector', 'regulatory_effect',
                       'regulatory_interaction_hash'])

    @staticmethod
    def build_network(genetic_network: pd.DataFrame,
                      regulatory_interaction: pd.DataFrame,
                      tf_gene: pd.DataFrame,
                      srna_interaction: pd.DataFrame,
                      tf: pd.DataFrame,
                      sigma: pd.DataFrame,
                      srna: pd.DataFrame,
                      gene: pd.DataFrame,
                      site: pd.DataFrame,
                      effector: pd.DataFrame,
                      conformation: pd.DataFrame,
                      conformation_effector: pd.DataFrame) -> pd.DataFrame:
        # ------------------
        # REGULATOR-GENE-TFBS
        # ------------------
        # GENETIC NETWORK
        # genetic_network drop regulated type == operon
        mask = genetic_network['regulated_type'] != 'operon'
        genetic_network = genetic_network[mask].reset_index(drop=True)

        # genetic_network_tf_* (regulator_id) + regulator_tf (tf_id)
        tf_mask = genetic_network['regulator_type'] == 'tf'
        genetic_network_tf = genetic_network[tf_mask].reset_index(drop=True)
        genetic_network_tf = pd.merge(genetic_network_tf, tf, left_on='regulator_id',
                                      right_on='transcription_factor_id')

        # genetic_network_sigma_* (regulator_id) + regulator_sigma (sigma_id)
        sigma_mask = genetic_network['regulator_type'] == 'sigma'
        genetic_network_sigma = genetic_network[sigma_mask].reset_index(drop=True)
        genetic_network_sigma = pd.merge(genetic_network_sigma, sigma, left_on='regulator_id', right_on='sigma_id')

        # vertical stacking of the resulting genetic_network
        genetic_network_tf = genetic_network_tf.reset_index(drop=True)
        genetic_network_sigma = genetic_network_sigma.reset_index(drop=True)
        genetic_network = pd.concat([genetic_network_tf, genetic_network_sigma])
        genetic_network = genetic_network.reset_index(drop=True)

        # genetic_network_*_gene (regulated_id) + gene (gene_id)
        gene_mask = genetic_network['regulated_type'] == 'gene'
        genetic_network_gene = genetic_network[gene_mask].reset_index(drop=True)
        genetic_network_gene = pd.merge(genetic_network_gene, gene, left_on='regulated_id', right_on='gene_id')

        # genetic_network_*_tf (regulated_name) + gene (gene_name)
        tf_mask = genetic_network['regulated_type'] == 'tf'
        genetic_network_tf = genetic_network[tf_mask].reset_index(drop=True)
        genetic_network_tf = pd.merge(genetic_network_tf, gene, left_on='regulated_name', right_on='gene_name')

        # genetic_network_*_sigma (regulated_id) + gene (gene_id)
        sigma_mask = genetic_network['regulated_type'] == 'sigma'
        genetic_network_sigma = genetic_network[sigma_mask].reset_index(drop=True)
        genetic_network_sigma = pd.merge(genetic_network_sigma, gene, left_on='regulated_id', right_on='gene_id')

        genetic_network_gene = genetic_network_gene.reset_index(drop=True)
        genetic_network_tf = genetic_network_tf.reset_index(drop=True)
        genetic_network_sigma = genetic_network_sigma.reset_index(drop=True)

        genetic_network = pd.concat([genetic_network_gene, genetic_network_tf, genetic_network_sigma])
        genetic_network = genetic_network.reset_index(drop=True)

        # reordering columns
        genetic_network = genetic_network.drop(columns=['regulated_id', 'regulated_name',
                                                        'regulator_type', 'regulated_type'])
        genetic_network = genetic_network.rename(columns={'function_interaction': 'regulatory_effect'})

        # REGULATORY INTERACTION
        # regulatory_interaction (conformation_id) + conformation (conformation_id)
        regulatory_interaction = pd.merge(regulatory_interaction, conformation, on='conformation_id')

        # regulatory_interaction (tf_id) + tf (tf_id)
        regulatory_interaction = pd.merge(regulatory_interaction, tf, on='transcription_factor_id')

        # regulatory_interaction (ri_first_gene_id) + gene (gene_id)
        regulatory_interaction = pd.merge(regulatory_interaction, gene, left_on='ri_first_gene_id', right_on='gene_id')

        # regulatory_interaction (site_id) + site (site_id): left
        regulatory_interaction = pd.merge(regulatory_interaction, site, how='left', on='site_id')

        # reordering columns
        regulatory_interaction = regulatory_interaction.assign(
            regulator_id=regulatory_interaction['transcription_factor_id'].copy()
        )
        regulatory_interaction = regulatory_interaction.drop(columns=['ri_first_gene_id'])
        regulatory_interaction = regulatory_interaction.rename(columns={'ri_function': 'regulatory_effect'})

        # TF-GENE
        # tf_gene (conformation_id) + conformation (conformation_id)
        tf_gene = pd.merge(tf_gene, conformation, on='conformation_id')

        # tf_gene (tf_id) + tf (tf_id)
        tf_gene = pd.merge(tf_gene, tf, on='transcription_factor_id')

        # tf_gene (object_id) + gene (gene_id)
        tf_gene = pd.merge(tf_gene, gene, left_on='object_id', right_on='gene_id')

        # tf_gene (site_id) + site (site_id): left
        tf_gene = pd.merge(tf_gene, site, how='left', on='site_id')

        # reordering columns
        tf_gene = tf_gene.assign(regulator_id=tf_gene['transcription_factor_id'].copy())
        tf_gene = tf_gene.drop(columns=['object_id'])
        tf_gene = tf_gene.rename(columns={'ri_function': 'regulatory_effect'})

        # SRNA
        # srna_interaction (srna_gene_id) + regulator (srna_gene_id)
        srna_interaction = pd.merge(srna_interaction, srna, on='srna_gene_id')

        # srna_interaction (srna_gene_regulated_id) + gene (gene_id)
        srna_interaction = pd.merge(srna_interaction, gene, left_on='srna_gene_regulated_id', right_on='gene_id')

        # reordering columns
        srna_interaction = srna_interaction.assign(regulator_id=srna_interaction['srna_gene_id'].copy())
        srna_interaction = srna_interaction.drop(columns=['srna_gene_regulated_id'])
        srna_interaction = srna_interaction.rename(columns={'srna_function': 'regulatory_effect'})

        # NETWORK
        # vertical stacking of all networks
        genetic_network = genetic_network.reset_index(drop=True)
        regulatory_interaction = regulatory_interaction.reset_index(drop=True)
        tf_gene = tf_gene.reset_index(drop=True)
        srna_interaction = srna_interaction.reset_index(drop=True)
        network = pd.concat([genetic_network, regulatory_interaction, tf_gene, srna_interaction])
        network = network.reset_index(drop=True)

        # ------------------
        # EFFECTOR-REGULATOR-GENE-TFBS
        # ------------------
        # EFFECTOR
        # network (conformation_id) + conformation_effector (conformation_id): left
        network = pd.merge(network, conformation_effector, how='left', on='conformation_id')
        network = pd.merge(network, effector, how='left', on='effector_id')
        return network

    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        network = network.assign(ncbi_taxonomy='511145')
        return network

    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        organism = select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)
        organism = organism.rename(columns={'protrend_id': 'organism'})
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = select_columns(regulator, 'protrend_id', 'regulator_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = select_columns(gene, 'protrend_id', 'gene_id')
        gene = gene.rename(columns={'protrend_id': 'gene'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = select_columns(tfbs, 'protrend_id', 'site_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform_effector(self, effector: pd.DataFrame) -> pd.DataFrame:
        effector = select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'effector'})
        return effector

    def regulondb_read(self, key: str):
        file, skiprows, columns = REGULONDB_STACK[key]
        reader = regulondb_reader(skiprows=skiprows, names=columns)
        return read(source=self.source, version=self.version, file=file, reader=reader,
                    default=pd.DataFrame(columns=columns))

    def transform(self) -> pd.DataFrame:
        tf = self.regulondb_read(key='tf')
        tf = transform_regulon_db_dataset(tf,
                                          selection=['transcription_factor_id'],
                                          duplicates=['transcription_factor_id'],
                                          nan=['transcription_factor_id'])

        sigma = self.regulondb_read(key='sigma')
        sigma = transform_regulon_db_dataset(sigma,
                                             selection=['sigma_id'],
                                             duplicates=['sigma_id'],
                                             nan=['sigma_id'])

        srna = self.regulondb_read(key='srna')
        srna = transform_regulon_db_dataset(srna,
                                            selection=['srna_gene_id'],
                                            duplicates=['srna_gene_id'],
                                            nan=['srna_gene_id'])

        target_gene = self.regulondb_read(key='target_gene')
        target_gene = transform_regulon_db_dataset(target_gene,
                                                   selection=['gene_id', 'gene_name'],
                                                   duplicates=['gene_id'],
                                                   nan=['gene_id'])

        site = self.regulondb_read(key='site')
        site = transform_regulon_db_dataset(site,
                                            selection=['site_id'],
                                            duplicates=['site_id'],
                                            nan=['site_id'])

        regulator_effector = self.regulondb_read(key='regulator_effector')
        regulator_effector = transform_regulon_db_dataset(regulator_effector,
                                                          selection=['effector_id'],
                                                          duplicates=['effector_id'],
                                                          nan=['effector_id'])

        conformation = self.regulondb_read(key='conformation')
        conformation = transform_regulon_db_dataset(conformation,
                                                    selection=['conformation_id', 'transcription_factor_id'],
                                                    duplicates=['conformation_id', 'transcription_factor_id'],
                                                    nan=['conformation_id', 'transcription_factor_id'])

        conformation_effector = self.regulondb_read(key='conformation_effector')
        conformation_effector = transform_regulon_db_dataset(conformation_effector,
                                                             selection=['conformation_id', 'effector_id'],
                                                             duplicates=['conformation_id', 'effector_id'],
                                                             nan=['conformation_id', 'effector_id'])

        genetic_network = self.regulondb_read(key='genetic_network')
        genetic_network = transform_regulon_db_dataset(genetic_network,
                                                       selection=['regulator_id', 'regulator_name', 'regulated_id',
                                                                  'regulated_name', 'function_interaction',
                                                                  'evidence', 'regulator_type', 'regulated_type'],
                                                       duplicates=['regulator_id', 'regulated_id',
                                                                   'function_interaction'],
                                                       nan=['regulator_id', 'regulated_id', 'function_interaction'])

        regulatory_interaction = self.regulondb_read(key='regulatory_interaction')
        regulatory_interaction = transform_regulon_db_dataset(regulatory_interaction,
                                                              selection=['conformation_id', 'site_id', 'ri_function',
                                                                         'ri_first_gene_id'],
                                                              duplicates=['conformation_id', 'site_id', 'ri_function',
                                                                          'ri_first_gene_id'],
                                                              nan=['conformation_id', 'ri_function',
                                                                   'ri_first_gene_id'])

        srna_interaction = self.regulondb_read(key='srna_interaction')
        srna_interaction = transform_regulon_db_dataset(srna_interaction,
                                                        selection=['srna_gene_id', 'srna_gene_regulated_id',
                                                                   'srna_function'],
                                                        duplicates=['srna_gene_id', 'srna_gene_regulated_id',
                                                                    'srna_function'],
                                                        nan=['srna_gene_id', 'srna_gene_regulated_id',
                                                             'srna_function'])

        tf_gene = self.regulondb_read(key='tf_gene')
        tf_gene = transform_regulon_db_dataset(tf_gene,
                                               selection=['conformation_id', 'object_id', 'site_id', 'ri_function'],
                                               duplicates=['conformation_id', 'object_id', 'site_id', 'ri_function'],
                                               nan=['conformation_id', 'object_id', 'ri_function'])

        network = self.build_network(genetic_network=genetic_network,
                                     regulatory_interaction=regulatory_interaction,
                                     tf_gene=tf_gene,
                                     srna_interaction=srna_interaction,
                                     tf=tf,
                                     sigma=sigma,
                                     srna=srna,
                                     gene=target_gene,
                                     site=site,
                                     effector=regulator_effector,
                                     conformation=conformation,
                                     conformation_effector=conformation_effector)

        organism = read_organism(source=self.source, version=self.version, columns=OrganismTransformer.columns)
        regulator = read_regulator(source=self.source, version=self.version, columns=RegulatorTransformer.columns)
        gene = read_gene(source=self.source, version=self.version, columns=GeneTransformer.columns)
        tfbs = read_tfbs(source=self.source, version=self.version, columns=TFBSTransformer.columns)
        effector = read_effector(source=self.source, version=self.version, columns=EffectorTransformer.columns)

        df = self._transform(network=network,
                             organism=organism, organism_key='ncbi_taxonomy',
                             regulator=regulator, regulator_key='regulator_id',
                             gene=gene, gene_key='gene_id',
                             tfbs=tfbs, tfbs_key='site_id',
                             effector=effector, effector_key='effector_id',
                             regulatory_effect_processor=regulatory_effect_regulondb)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToEffectorConnector(RegulonDBConnector,
                                               source='regulondb',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Effector,
                                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='effector')
        self.stack_connections(df)


class RegulatoryInteractionToRegulatorConnector(RegulonDBConnector,
                                                source='regulondb',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_connections(df)


class RegulatoryInteractionToGeneConnector(RegulonDBConnector,
                                           source='regulondb',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_connections(df)


class RegulatoryInteractionToTFBSConnector(RegulonDBConnector,
                                           source='regulondb',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='tfbs')
        self.stack_connections(df)


class RegulatorToEffectorConnector(RegulonDBConnector,
                                   source='regulondb',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Effector,
                                   register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='effector')
        self.stack_connections(df)


class RegulatorToGeneConnector(RegulonDBConnector,
                               source='regulondb',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_connections(df)


class RegulatorToTFBSConnector(RegulonDBConnector,
                               source='regulondb',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='tfbs')
        self.stack_connections(df)


class GeneToTFBSConnector(RegulonDBConnector,
                          source='regulondb',
                          version='0.0.0',
                          from_node=Gene,
                          to_node=TFBS,
                          register=True):

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_connections(df)
