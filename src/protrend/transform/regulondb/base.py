from typing import List

import pandas as pd

from protrend.io import read_from_stack, read_txt
from protrend.transform import Transformer, Connector
from protrend.transform.processors import to_set_list, apply_processors, to_list
from protrend.utils import SetList


class RegulondbTransformer(Transformer):
    default_source = 'regulondb'
    default_version = '0.0.0'

    regulondb_stack = {'site': 'site.txt',
                       'ri': 'regulatory_interaction.txt',
                       'tf_gene_interaction': 'tf_gene_interaction.txt',
                       'gene': 'gene.txt',
                       'tf': 'transcription_factor.txt',
                       'srna': 'srna_interaction.txt',
                       'sigma': 'sigma_tmp.txt',
                       'tu': 'transcription_unit.txt'}

    tfbs_columns = SetList(['site_id', 'site_posleft', 'site_posright', 'site_sequence', 'site_note',
                            'site_internal_comment', 'key_id_org', 'site_length'])

    tf_columns = SetList(['transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                          'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                          'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence'])
    srna_columns = SetList(['srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id',
                            'srna_function', 'srna_posleft', 'srna_posright', 'srna_sequence',
                            'srna_regulation_type', 'srna_mechanis', 'srna_note'])
    sigma_columns = SetList(['sigma_id', 'sigma_name', 'sigma_synonyms', 'sigma_gene_id', 'sigma_gene_name',
                             'sigma_coregulators', 'sigma_notes', 'sigma_sigmulon_genes', 'key_id_org'])

    effector_columns = SetList(['effector_id', 'effector_name', 'category', 'effector_type', 'effector_note',
                                'effector_internal_comment', 'key_id_org'])

    conformation_columns = SetList(['conformation_id', 'transcription_factor_id', 'final_state', 'conformation_note',
                                    'interaction_type', 'conformation_internal_comment', 'key_id_org',
                                    'conformation_type', 'apo_holo_conformation'])

    promoter_columns = SetList(['promoter_id', 'promoter_name', 'promoter_strand', 'pos_1', 'sigma_factor',
                                'basal_trans_val', 'equilibrium_const', 'kinetic_const', 'strength_seq',
                                'promoter_sequence', 'key_id_org', 'promoter_note', 'promoter_internal_comment'])

    conformation_eff_columns = SetList(['conformation_id', 'effector_id'])

    gene_columns = SetList(['gene_id', 'gene_name', 'gene_posleft', 'gene_posright', 'gene_strand',
                            'gene_sequence', 'gc_content', 'cri_score', 'gene_note',
                            'gene_internal_comment', 'key_id_org', 'gene_type'])

    tu_columns = SetList(['transcription_unit_id', 'promoter_id', 'transcription_unit_name', 'operon_id',
                          'key_id_org', 'transcription_unit_note', 'tu_internal_comment'])

    tu_gene_columns = SetList(['transcription_unit_id', 'gene_id'])

    gen_net_columns = SetList(['regulator_id', 'regulator_name', 'regulated_id', 'regulated_name',
                               'function_interaction', 'evidence', 'regulator_type', 'regulated_type'])

    ri_columns = SetList(['regulatory_interaction_id', 'conformation_id', 'promoter_id', 'site_id',
                          'ri_function', 'center_position', 'ri_dist_first_gene', 'ri_first_gene_id',
                          'affinity_exp', 'regulatory_interaction_note', 'ri_internal_comment', 'key_id_org',
                          'ri_sequence', 'ri_orientation', 'ri_sequence_orientation'])

    tf_gene_columns = SetList(['regulatory_interaction_id', 'conformation_id', 'object_id', 'site_id',
                               'ri_function', 'center_position', 'ri_dist_first_gene', 'ri_first_gene_id',
                               'affinity_exp', 'regulatory_interaction_note', 'ri_internal_comment', 'key_id_org',
                               'ri_sequence', 'ri_orientation', 'ri_sequence_orientation', 'object_type'])

    def _build(self, df: pd.DataFrame, selection: List[str], duplicates: List[str], nan: List[str]) -> pd.DataFrame:
        df = self.select_columns(df, *selection)
        df = self.drop_duplicates(df, subset=duplicates, perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=nan)
        return df

    def _build_tfbs(self):
        df = read_from_stack(stack=self.regulondb_stack, file='site', default_columns=self.tfbs_columns,
                             reader=read_txt, skiprows=35, names=self.tfbs_columns)
        return self._build(df=df, selection=['site_id'], duplicates=['site_id'], nan=['site_id'])

    def _build_gene(self):
        df = read_from_stack(stack=self.regulondb_stack, file='gene',
                             default_columns=self.gene_columns, reader=read_txt,
                             skiprows=39, names=self.gene_columns)
        return self._build(df=df, selection=['gene_id', 'gene_name'], duplicates=['gene_id'], nan=['gene_id'])

    def _build_tu(self):
        df = read_from_stack(stack=self.regulondb_stack, file='tu',
                             default_columns=self.tu_columns, reader=read_txt,
                             skiprows=34, names=self.tu_columns)
        return self._build(df=df,
                           selection=['transcription_unit_id', 'promoter_id', 'transcription_unit_name', 'operon_id'],
                           duplicates=['transcription_unit_id'],
                           nan=['transcription_unit_id'])

    def _build_operon(self):
        tu = self._build_tu()

        tu_gene = read_from_stack(stack=self.regulondb_stack, file='tu_gene',
                                  default_columns=self.tu_gene_columns, reader=read_txt,
                                  skiprows=29, names=self.tu_gene_columns)
        df = pd.merge(tu, tu_gene, on='transcription_unit_id')

        return self.group_by(df, column='operon_id', aggregation={}, default=to_set_list)

    def _build_tf(self):
        df = read_from_stack(stack=self.regulondb_stack, file='tf',
                             default_columns=self.tf_columns, reader=read_txt,
                             skiprows=38, names=self.tf_columns)
        return self._build(df=df,
                           selection=['transcription_factor_id', 'transcription_factor_name'],
                           duplicates=['transcription_factor_id'],
                           nan=['transcription_factor_id'])

    def _build_srna(self):
        df = read_from_stack(stack=self.regulondb_stack, file='srna',
                             default_columns=self.srna_columns, reader=read_txt,
                             skiprows=38, names=self.srna_columns)
        return self._build(df=df,
                           selection=['srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id'],
                           duplicates=['srna_gene_id', 'srna_gene_regulated_id'],
                           nan=['srna_gene_id', 'srna_gene_regulated_id'])

    def _build_sigma(self):
        df = read_from_stack(stack=self.regulondb_stack, file='sigma',
                             default_columns=self.sigma_columns, reader=read_txt,
                             skiprows=36, names=self.sigma_columns)
        return self._build(df=df,
                           selection=['sigma_id', 'sigma_name', 'sigma_gene_id', 'sigma_gene_name'],
                           duplicates=['sigma_id'],
                           nan=['sigma_id'])

    def _build_effector(self):
        df = read_from_stack(stack=self.regulondb_stack, file='effector', default_columns=self.effector_columns,
                             reader=read_txt, skiprows=34, names=self.effector_columns)
        return self._build(df=df,
                           selection=['effector_id', 'effector_name'],
                           duplicates=['effector_id'],
                           nan=['effector_id'])

    def _build_conformation(self):
        df = read_from_stack(stack=self.regulondb_stack, file='conformation',
                             default_columns=self.conformation_columns,
                             reader=read_txt, skiprows=36, names=self.conformation_columns)

        return self._build(df=df,
                           selection=['conformation_id', 'transcription_factor_id'],
                           duplicates=['conformation_id'],
                           nan=['conformation_id'])

    def _build_promoter(self):
        df = read_from_stack(stack=self.regulondb_stack, file='promoter',
                             default_columns=self.promoter_columns, reader=read_txt,
                             skiprows=40, names=self.promoter_columns)
        return self._build(df=df,
                           selection=['promoter_id', 'promoter_name', 'promoter_sequence'],
                           duplicates=['promoter_id'],
                           nan=['promoter_id', 'promoter_sequence'])

    def _build_conformation_effector(self):
        df = read_from_stack(stack=self.regulondb_stack, file='conformation_effector',
                             default_columns=self.conformation_eff_columns,
                             reader=read_txt, skiprows=29, names=self.conformation_eff_columns)
        return self._build(df=df,
                           selection=['conformation_id', 'effector_id'],
                           duplicates=['conformation_id', 'effector_id'],
                           nan=['conformation_id', 'effector_id'])

    def _build_regulator_effector(self):
        tf = self._build_tf()
        conformation = self._build_conformation()

        tf_conformation = pd.merge(tf, conformation, on='transcription_factor_id')

        effector = self._build_effector()
        conformation_effector = self._build_conformation_effector()

        conformation_effector = pd.merge(effector, conformation_effector, on='effector_id')

        tf_effector = pd.merge(tf_conformation, conformation_effector, on='conformation_id')
        return tf_effector

    def _build_genetic_network(self):
        df = read_from_stack(stack=self.regulondb_stack, file='genetic_network',
                             default_columns=self.gen_net_columns, reader=read_txt,
                             skiprows=35, names=self.gen_net_columns)
        return self._build(df=df,
                           selection=['regulator_id', 'regulator_name', 'regulated_id' 'regulated_name',
                                      'function_interaction', 'regulator_type', 'regulated_type'],
                           duplicates=['regulator_id', 'regulated_id'],
                           nan=['regulator_id', 'regulated_id'])

    def _build_regulatory_interaction(self):
        df = read_from_stack(stack=self.regulondb_stack, file='ri', default_columns=self.ri_columns,
                             reader=read_txt, skiprows=42, names=self.ri_columns)
        return self._build(df=df,
                           selection=['regulatory_interaction_id', 'conformation_id', 'promoter_id', 'site_id',
                                      'ri_dist_first_gene'],
                           duplicates=['conformation_id', 'site_id', 'ri_dist_first_gene'],
                           nan=['conformation_id', 'ri_dist_first_gene'])

    def _build_tf_gene_interaction(self):
        df = read_from_stack(stack=self.regulondb_stack, file='tf_gene_interaction',
                             default_columns=self.tf_gene_columns, reader=read_txt,
                             skiprows=43, names=self.tf_gene_columns)
        return self._build(df=df,
                           selection=['regulatory_interaction_id', 'conformation_id', 'object_id',
                                      'site_id', 'ri_first_gene_id'],
                           duplicates=['conformation_id', 'object_id', 'site_id'],
                           nan=['conformation_id', 'object_id', 'ri_first_gene_id'])

    def _build_regulator_gene(self):
        tf = self._build_tf()
        conformation = self._build_conformation()
        tf_conformation = pd.merge(conformation, tf, on='transcription_factor_id')

        sigma = self._build_sigma()
        gene = self._build_gene()

        # Genetic network
        gen_net = self._build_genetic_network()

        tf_gen_net = pd.merge(gen_net, tf, left_on='regulator_id', right_on='transcription_factor_id')
        sigma_gen_net = pd.merge(gen_net, sigma, left_on='regulator_id', right_on='sigma_id')

        reg_gen_net = pd.concat([tf_gen_net, sigma_gen_net], axis=0)
        reg_gen_net = pd.merge(reg_gen_net, gene, left_on='regulated_id', right_on='gene_id')
        reg_gen_net = self.select_columns('regulator_id', 'regulated_id')
        reg_gen_net = reg_gen_net.rename(columns={'regulated_id': 'gene_id'})

        # Regulatory Interaction
        ri = self._build_regulatory_interaction()
        tf_ri = pd.merge(ri, tf_conformation, on='conformation_id')
        tf_gene_ri = pd.merge(tf_ri, gene, left_on='ri_dist_first_gene', right_on='gene_id')
        tf_gene_ri = self.select_columns('transcription_factor_id', 'gene_id')
        tf_gene_ri = tf_gene_ri.rename(columns={'transcription_factor_id': 'regulator_id'})

        # TF-Gene Interaction
        tf_gene_i = self._build_tf_gene_interaction()
        tf_tf_gene_i = pd.merge(tf_gene_i, tf_conformation, on='conformation_id')
        tf_gene_tf_gene_i = pd.merge(tf_tf_gene_i, gene, left_on='object_id', right_on='gene_id')
        tf_gene_tf_gene_i = self.select_columns('transcription_factor_id', 'gene_id')
        tf_gene_tf_gene_i = tf_gene_tf_gene_i.rename(columns={'transcription_factor_id': 'regulator_id'})

        # SRNA Interaction
        srna = self._build_srna()
        srna = self.select_columns('srna_gene_id', 'srna_gene_regulated_id')
        srna = srna.rename(columns={'srna_gene_id': 'regulator_id', 'srna_gene_regulated_id': 'gene_id'})

        reg_gen = pd.concat([reg_gen_net, tf_gene_ri, tf_gene_tf_gene_i, srna], axis=0)
        reg_gen = self.drop_duplicates(reg_gen, subset=['regulator_id', 'gene_id'],
                                       perfect_match=True, preserve_nan=True)
        reg_gen = reg_gen.dropna(subset=['regulator_id', 'gene_id'])

        return reg_gen

    def _build_regulator_operon(self):
        tf = self._build_tf()
        sigma = self._build_sigma()
        operon = self._build_operon()

        # Genetic network
        gen_net = self._build_genetic_network()

        tf_gen_net = pd.merge(gen_net, tf, left_on='regulator_id', right_on='transcription_factor_id')
        sigma_gen_net = pd.merge(gen_net, sigma, left_on='regulator_id', right_on='sigma_id')

        reg_op_net = pd.concat([tf_gen_net, sigma_gen_net], axis=0)
        reg_op_net = pd.merge(reg_op_net, operon, left_on='regulated_id', right_on='operon_id')
        reg_op_net = self.select_columns('regulator_id', 'regulated_id')
        reg_op_net = reg_op_net.rename(columns={'regulated_id': 'operon_id'})

        # All Regulator-Gene interactions
        reg_gen = self._build_regulator_gene()
        op_by_gene = apply_processors(operon, gene_id=to_list)
        op_by_gene = operon.explode(column='gene_id')
        reg_gen_op = pd.merge(reg_gen, op_by_gene, on='gene_id')
        reg_gen_op = self.select_columns('regulator_id', 'operon_id')

        reg_op = pd.concat([reg_op_net, reg_gen_op], axis=0)
        reg_op = self.drop_duplicates(reg_op, subset=['regulator_id', 'operon_id'],
                                      perfect_match=True, preserve_nan=True)
        reg_op = reg_op.dropna(subset=['regulator_id', 'operon_id'])

        return reg_op

    def _build_operon_gene(self):
        operon = self._build_operon()
        operon_gene = apply_processors(operon, gene_id=to_list)
        operon_gene = operon.explode(column='gene_id')
        operon_gene = self.select_columns('operon_id', 'gene_id')

        return operon_gene

    def _build_operon_promoter(self):
        operon = self._build_operon()
        operon_promoter = apply_processors(operon, promoter_id=to_list)
        operon_promoter = operon.explode(column='promoter_id')
        operon_promoter = self.select_columns('operon_id', 'promoter_id')

        return operon_promoter

    def _build_operon_tfbs(self):
        tf = self._build_tf()
        tfbs = self._build_tfbs()
        conformation = self._build_conformation()
        tf_conformation = pd.merge(conformation, tf, on='transcription_factor_id')

        # Regulatory Interaction
        ri = self._build_regulatory_interaction()
        tf_ri = pd.merge(ri, tf_conformation, on='conformation_id')
        tf_ri = self.select_columns('transcription_factor_id', 'site_id')
        tf_ri = self.drop_duplicates(tf_ri, subset=['transcription_factor_id', 'site_id'],
                                     perfect_match=True, preserve_nan=True)
        tf_ri = tf_ri.dropna(subset=['transcription_factor_id', 'site_id'])
        tf_ri = tf_ri.rename(columns={'transcription_factor_id': 'regulator_id'})

        # TF-Gene Interaction
        tf_gene_i = self._build_tf_gene_interaction()
        tf_tf_gene_i = pd.merge(tf_gene_i, tf_conformation, on='conformation_id')
        tf_gene_i = self.select_columns('transcription_factor_id', 'site_id')
        tf_gene_i = self.drop_duplicates(tf_gene_i, subset=['transcription_factor_id', 'site_id'],
                                         perfect_match=True, preserve_nan=True)
        tf_gene_i = tf_gene_i.dropna(subset=['transcription_factor_id', 'site_id'])
        tf_gene_i = tf_gene_i.rename(columns={'transcription_factor_id': 'regulator_id'})

        tf_tfbs = pd.concat([tf_ri, tf_gene_i], axis=0)
        tf_tfbs = pd.merge(tf_tfbs, tfbs, on='site_id')

        tf_operon = self._build_regulator_operon()

        operon_tfbs = pd.merge(tf_operon, tf_tfbs, on='regulator_id')
        operon_tfbs = self.select_columns('operon_id', 'site_id')
        operon_tfbs = self.drop_duplicates(operon_tfbs, subset=['operon_id', 'site_id'],
                                           perfect_match=True, preserve_nan=True)
        operon_tfbs = operon_tfbs.dropna(subset=['operon_id', 'site_id'])

        return operon_tfbs

    def _build_gene_tfbs(self):
        operon = self._build_operon()
        operon_tfbs = self._build_operon_tfbs()

        gene_tfbs = pd.merge(operon_tfbs, operon, on='operon_id')
        gene_tfbs = apply_processors(gene_tfbs, gene_id=to_list)
        gene_tfbs = gene_tfbs.explode(column='gene_id')

        gene_tfbs = self.select_columns('gene_id', 'site_id')
        gene_tfbs = self.drop_duplicates(gene_tfbs, subset=['gene_id', 'site_id'],
                                         perfect_match=True, preserve_nan=True)
        gene_tfbs = gene_tfbs.dropna(subset=['gene_id', 'site_id'])

        return gene_tfbs


class RegulondbConnector(Connector):
    default_source: str = 'regulondb'
    default_version: str = '0.0.0'
