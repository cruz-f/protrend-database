from typing import List

import pandas as pd

from protrend.io import read_json_frame, read_from_stack
from protrend.model import RegulatoryInteraction, Regulator, TFBS, Gene, Effector
from protrend.transform import RegulatoryInteractionMixIn, Transformer
from protrend.transform.regulondb.base import RegulondbTransformer, RegulondbConnector, regulondb_reader
from protrend.transform.regulondb.effector import EffectorTransformer
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.regulondb.organism import OrganismTransformer
from protrend.transform.regulondb.regulator import RegulatorTransformer
from protrend.transform.regulondb.tfbs import TFBSTransformer
from protrend.utils import SetList
from protrend.utils.processors import (apply_processors, regulatory_effect_regulondb, to_int_str)


def transform_regulon_db_dataset(df: pd.DataFrame,
                                 selection: List[str],
                                 duplicates: List[str],
                                 nan: List[str]) -> pd.DataFrame:
    df = Transformer.select_columns(df, *selection)
    df = df.dropna(subset=nan)
    df = Transformer.drop_empty_string(df, *nan)
    df = Transformer.drop_duplicates(df, subset=duplicates, perfect_match=True)
    return df


class RegulatoryInteractionTransformer(RegulatoryInteractionMixIn, RegulondbTransformer,
                                       source='regulondb',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=80,
                                       register=True):
    default_transform_stack = {'organism': 'integrated_organism.json',
                               'effector': 'integrated_effector.json',
                               'regulator': 'integrated_regulator.json',
                               'gene': 'integrated_gene.json',
                               'tfbs': 'integrated_tfbs.json',
                               'genetic_network': 'genetic_network.txt',
                               'regulatory_interaction': 'regulatory_interaction.txt',
                               'srna_interaction': 'srna_interaction.txt',
                               'tf_gene': 'tf_gene_interaction.txt',
                               'tf': 'transcription_factor.txt',
                               'sigma': 'sigma_tmp.txt',
                               'srna': 'srna_interaction.txt',
                               'target_gene': 'gene.txt',
                               'site': 'site.txt',
                               'regulator_effector': 'effector.txt',
                               'conformation': 'conformation.txt',
                               'conformation_effector': 'conformation_effector_link.txt'}

    genetic_network_columns = SetList(['regulator_id', 'regulator_name', 'regulated_id', 'regulated_name',
                                       'function_interaction', 'evidence', 'regulator_type', 'regulated_type'])
    regulatory_interaction_columns = SetList(['regulatory_interaction_id', 'conformation_id', 'promoter_id', 'site_id',
                                              'ri_function', 'center_position', 'ri_dist_first_gene',
                                              'ri_first_gene_id', 'affinity_exp', 'regulatory_interaction_note',
                                              'ri_internal_comment', 'key_id_org', 'ri_sequence', 'ri_orientation',
                                              'ri_sequence_orientation'])
    srna_interaction_columns = SetList(['srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id',
                                        'srna_function', 'srna_posleft', 'srna_posright', 'srna_sequence',
                                        'srna_regulation_type', 'srna_mechanis', 'srna_note'])
    tf_gene_columns = SetList(['regulatory_interaction_id', 'conformation_id', 'object_id', 'site_id',
                               'ri_function', 'center_position', 'ri_dist_first_gene', 'ri_first_gene_id',
                               'affinity_exp', 'regulatory_interaction_note', 'ri_internal_comment', 'key_id_org',
                               'ri_sequence', 'ri_orientation', 'ri_sequence_orientation', 'object_type'])
    tf_columns = SetList(['transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                          'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                          'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence'])
    sigma_columns = SetList(['sigma_id', 'sigma_name', 'sigma_synonyms', 'sigma_gene_id', 'sigma_gene_name',
                             'sigma_coregulators', 'sigma_notes', 'sigma_sigmulon_genes', 'key_id_org'])
    srna_columns = SetList(['srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id',
                            'srna_function', 'srna_posleft', 'srna_posright', 'srna_sequence',
                            'srna_regulation_type', 'srna_mechanis', 'srna_note'])
    target_gene_columns = SetList(['gene_id', 'gene_name', 'gene_posleft', 'gene_posright', 'gene_strand',
                                   'gene_sequence', 'gc_content', 'cri_score', 'gene_note',
                                   'gene_internal_comment', 'key_id_org', 'gene_type'])
    site_columns = SetList(['site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
                            'site_internal_comment', 'key_id_org', 'site_length'])
    regulator_effector_columns = SetList(['effector_id', 'effector_name', 'category', 'effector_type', 'effector_note',
                                          'effector_internal_comment', 'key_id_org'])
    conformation_columns = SetList(['conformation_id', 'transcription_factor_id', 'final_state', 'conformation_note',
                                    'interaction_type', 'conformation_internal_comment', 'key_id_org',
                                    'conformation_type', 'apo_holo_conformation'])
    conformation_effector_columns = SetList(['effector_id', 'conformation_id'])

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
        organism = self.select_columns(organism, 'protrend_id', 'ncbi_taxonomy')
        organism = apply_processors(organism, ncbi_taxonomy=to_int_str)
        organism = organism.rename(columns={'protrend_id': 'organism'})
        return organism

    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        regulator = self.select_columns(regulator, 'protrend_id', 'regulator_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator'})
        return regulator

    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = self.select_columns(gene, 'protrend_id', 'gene_id')
        gene = gene.rename(columns={'protrend_id': 'gene'})
        return gene

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        tfbs = self.select_columns(tfbs, 'protrend_id', 'site_id')
        tfbs = tfbs.rename(columns={'protrend_id': 'tfbs'})
        return tfbs

    def transform_effector(self, effector: pd.DataFrame) -> pd.DataFrame:
        effector = self.select_columns(effector, 'protrend_id', 'effector_id')
        effector = effector.rename(columns={'protrend_id': 'effector'})
        return effector

    def transform(self) -> pd.DataFrame:
        tf_reader = regulondb_reader(skiprows=38, names=self.tf_columns)
        tf = read_from_stack(stack=self.transform_stack, key='tf',
                             columns=self.tf_columns, reader=tf_reader)
        tf = transform_regulon_db_dataset(tf,
                                          selection=['transcription_factor_id'],
                                          duplicates=['transcription_factor_id'],
                                          nan=['transcription_factor_id'])

        sigma_reader = regulondb_reader(skiprows=36, names=self.sigma_columns)
        sigma = read_from_stack(stack=self.transform_stack, key='sigma',
                                columns=self.sigma_columns, reader=sigma_reader)
        sigma = transform_regulon_db_dataset(sigma,
                                             selection=['sigma_id'],
                                             duplicates=['sigma_id'],
                                             nan=['sigma_id'])

        srna_reader = regulondb_reader(skiprows=38, names=self.srna_columns)
        srna = read_from_stack(stack=self.transform_stack, key='srna',
                               columns=self.srna_columns, reader=srna_reader)
        srna = transform_regulon_db_dataset(srna,
                                            selection=['srna_gene_id'],
                                            duplicates=['srna_gene_id'],
                                            nan=['srna_gene_id'])

        target_gene_reader = regulondb_reader(skiprows=39, names=self.target_gene_columns)
        target_gene = read_from_stack(stack=self.transform_stack, key='target_gene',
                                      columns=self.target_gene_columns, reader=target_gene_reader)
        target_gene = transform_regulon_db_dataset(target_gene,
                                                   selection=['gene_id', 'gene_name'],
                                                   duplicates=['gene_id'],
                                                   nan=['gene_id'])

        site_reader = regulondb_reader(skiprows=35, names=self.site_columns)
        site = read_from_stack(stack=self.transform_stack, key='site',
                               columns=self.site_columns, reader=site_reader)
        site = transform_regulon_db_dataset(site,
                                            selection=['site_id'],
                                            duplicates=['site_id'],
                                            nan=['site_id'])

        effector_reader = regulondb_reader(skiprows=34, names=self.regulator_effector_columns)
        regulator_effector = read_from_stack(stack=self.transform_stack, key='regulator_effector',
                                             columns=self.regulator_effector_columns, reader=effector_reader)
        regulator_effector = transform_regulon_db_dataset(regulator_effector,
                                                          selection=['effector_id'],
                                                          duplicates=['effector_id'],
                                                          nan=['effector_id'])

        conformation_reader = regulondb_reader(skiprows=36, names=self.conformation_columns)
        conformation = read_from_stack(stack=self.transform_stack, key='conformation',
                                       columns=self.conformation_columns, reader=conformation_reader)
        conformation = transform_regulon_db_dataset(conformation,
                                                    selection=['conformation_id', 'transcription_factor_id'],
                                                    duplicates=['conformation_id', 'transcription_factor_id'],
                                                    nan=['conformation_id', 'transcription_factor_id'])

        conformation_effector_reader = regulondb_reader(skiprows=29,
                                                        names=self.conformation_effector_columns)
        conformation_effector = read_from_stack(stack=self.transform_stack, key='conformation_effector',
                                                columns=self.conformation_effector_columns,
                                                reader=conformation_effector_reader)
        conformation_effector = transform_regulon_db_dataset(conformation_effector,
                                                             selection=['conformation_id', 'effector_id'],
                                                             duplicates=['conformation_id', 'effector_id'],
                                                             nan=['conformation_id', 'effector_id'])

        genetic_network_reader = regulondb_reader(skiprows=43, names=self.genetic_network_columns)
        genetic_network = read_from_stack(stack=self.transform_stack, key='genetic_network',
                                          columns=self.genetic_network_columns,
                                          reader=genetic_network_reader)
        genetic_network = transform_regulon_db_dataset(genetic_network,
                                                       selection=['regulator_id', 'regulator_name', 'regulated_id',
                                                                  'regulated_name', 'function_interaction',
                                                                  'evidence', 'regulator_type', 'regulated_type'],
                                                       duplicates=['regulator_id', 'regulated_id',
                                                                   'function_interaction'],
                                                       nan=['regulator_id', 'regulated_id', 'function_interaction'])

        regulatory_interaction_reader = regulondb_reader(skiprows=42, names=self.regulatory_interaction_columns)
        regulatory_interaction = read_from_stack(stack=self.transform_stack, key='regulatory_interaction',
                                                 columns=self.regulatory_interaction_columns,
                                                 reader=regulatory_interaction_reader)
        regulatory_interaction = transform_regulon_db_dataset(regulatory_interaction,
                                                              selection=['conformation_id', 'site_id', 'ri_function',
                                                                         'ri_first_gene_id'],
                                                              duplicates=['conformation_id', 'site_id', 'ri_function',
                                                                          'ri_first_gene_id'],
                                                              nan=['conformation_id', 'ri_function',
                                                                   'ri_first_gene_id'])

        srna_interaction_reader = regulondb_reader(skiprows=38, names=self.srna_interaction_columns)
        srna_interaction = read_from_stack(stack=self.transform_stack, key='srna_interaction',
                                           columns=self.srna_interaction_columns,
                                           reader=srna_interaction_reader)
        srna_interaction = transform_regulon_db_dataset(srna_interaction,
                                                        selection=['srna_gene_id', 'srna_gene_regulated_id',
                                                                   'srna_function'],
                                                        duplicates=['srna_gene_id', 'srna_gene_regulated_id',
                                                                    'srna_function'],
                                                        nan=['srna_gene_id', 'srna_gene_regulated_id',
                                                             'srna_function'])

        tf_gene_reader = regulondb_reader(skiprows=43, names=self.tf_gene_columns)
        tf_gene = read_from_stack(stack=self.transform_stack, key='tf_gene',
                                  columns=self.tf_gene_columns,
                                  reader=tf_gene_reader)
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

        # noinspection DuplicatedCode
        organism = read_from_stack(stack=self.transform_stack, key='organism',
                                   columns=OrganismTransformer.columns, reader=read_json_frame)
        regulator = read_from_stack(stack=self.transform_stack, key='regulator',
                                    columns=RegulatorTransformer.columns, reader=read_json_frame)
        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=GeneTransformer.columns, reader=read_json_frame)
        tfbs = read_from_stack(stack=self.transform_stack, key='tfbs',
                               columns=TFBSTransformer.columns, reader=read_json_frame)
        effector = read_from_stack(stack=self.transform_stack, key='effector',
                                   columns=EffectorTransformer.columns, reader=read_json_frame)

        df = self._transform(network=network,
                             organism=organism, organism_key='ncbi_taxonomy',
                             regulator=regulator, regulator_key='regulator_id',
                             gene=gene, gene_key='gene_id',
                             tfbs=tfbs, tfbs_key='site_id',
                             effector=effector, effector_key='effector_id',
                             regulatory_effect_processor=regulatory_effect_regulondb)

        self.stack_transformed_nodes(df)
        return df


class RegulatoryInteractionToEffectorConnector(RegulondbConnector,
                                               source='regulondb',
                                               version='0.0.0',
                                               from_node=RegulatoryInteraction,
                                               to_node=Effector,
                                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='effector')
        self.stack_json(df)


class RegulatoryInteractionToRegulatorConnector(RegulondbConnector,
                                                source='regulondb',
                                                version='0.0.0',
                                                from_node=RegulatoryInteraction,
                                                to_node=Regulator,
                                                register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='regulator')
        self.stack_json(df)


class RegulatoryInteractionToGeneConnector(RegulondbConnector,
                                           source='regulondb',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=Gene,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='gene')
        self.stack_json(df)


class RegulatoryInteractionToTFBSConnector(RegulondbConnector,
                                           source='regulondb',
                                           version='0.0.0',
                                           from_node=RegulatoryInteraction,
                                           to_node=TFBS,
                                           register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    target_column='tfbs')
        self.stack_json(df)


class RegulatorToEffectorConnector(RegulondbConnector,
                                   source='regulondb',
                                   version='0.0.0',
                                   from_node=Regulator,
                                   to_node=Effector,
                                   register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='effector')
        self.stack_json(df)


class RegulatorToGeneConnector(RegulondbConnector,
                               source='regulondb',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=Gene,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='gene')
        self.stack_json(df)


class RegulatorToTFBSConnector(RegulondbConnector,
                               source='regulondb',
                               version='0.0.0',
                               from_node=Regulator,
                               to_node=TFBS,
                               register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='regulator', target_column='tfbs')
        self.stack_json(df)


class GeneToTFBSConnector(RegulondbConnector,
                          source='regulondb',
                          version='0.0.0',
                          from_node=Gene,
                          to_node=TFBS,
                          register=True):
    default_connect_stack = {'rin': 'integrated_regulatoryinteraction.json'}

    def connect(self):
        df = self.create_connection(source='rin', target='rin',
                                    source_column='gene', target_column='tfbs')
        self.stack_json(df)
