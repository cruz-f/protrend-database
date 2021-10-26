import pandas as pd

from protrend.io import read_from_stack, read_txt, read_json_frame
from protrend.model.model import Regulator
from protrend.transform.processors import apply_processors, rstrip, lstrip
from protrend.transform.regulondb.base import RegulondbTransformer
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.utils import SetList


class RegulatorTransformer(RegulondbTransformer):
    default_node = Regulator
    default_transform_stack = {'tf': 'transcription_factor.txt',
                               'srna': 'srna_interaction.txt',
                               'sigma': 'sigma_tmp.txt',
                               'gene': 'integrated_gene.json'}
    default_order = 90
    columns = SetList(['locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession', 'sequence',
                       'strand', 'start', 'stop', 'mechanism',
                       'transcription_factor_id', 'transcription_factor_name', 'transcription_factor_family',
                       'srna_id', 'srna_gene_id', 'sigma_id', 'sigma_name', 'sigma_id', 'sigma_name',
                       'gene_id', 'gene_name_lower', 'gene_protrend_id', 'protrend_id'])
    tf_columns = SetList(['transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                          'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                          'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence'])
    srna_columns = SetList(['srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id',
                            'srna_function', 'srna_posleft', 'srna_posright', 'srna_sequence',
                            'srna_regulation_type', 'srna_mechanis', 'srna_note'])
    sigma_columns = SetList(['sigma_id', 'sigma_name', 'sigma_synonyms', 'sigma_gene_id', 'sigma_gene_name',
                             'sigma_coregulators', 'sigma_notes', 'sigma_sigmulon_genes', 'key_id_org'])

    def _transform_tf(self, tf: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        tf = apply_processors(tf, transcription_factor_name=[rstrip, lstrip])
        tf = tf.dropna(subset=['transcription_factor_name'])
        tf = self.drop_duplicates(df=tf, subset=['transcription_factor_name'], perfect_match=True, preserve_nan=True)

        tf['mechanism'] = 'transcription factor'
        tf['gene_name_lower'] = tf['transcription_factor_name'].str.lower()

        tf = pd.merge(tf, gene, on='gene_name_lower')

        return tf

    def _transform_srna(self, srna: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        srna = apply_processors(srna, srna_gene_id=[rstrip, lstrip])
        srna = srna.dropna(subset=['srna_gene_id'])
        srna = self.drop_duplicates(df=srna, subset=['srna_gene_id'], perfect_match=True, preserve_nan=True)

        srna['mechanism'] = 'small RNA (sRNA)'

        srna = pd.merge(srna, gene, left_on='srna_gene_id', right_on='gene_id')

        return srna

    def _transform_sigma(self, sigma: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        sigma = apply_processors(sigma, sigma_gene_id=[rstrip, lstrip])
        sigma = sigma.dropna(subset=['sigma_gene_id'])
        sigma = self.drop_duplicates(df=sigma, subset=['sigma_gene_id'], perfect_match=True, preserve_nan=True)

        sigma['mechanism'] = 'sigma factor'

        sigma = pd.merge(sigma, gene, left_on='sigma_gene_id', right_on='gene_id')

        return sigma

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.tf_columns, reader=read_txt,
                             skiprows=38, names=self.tf_columns)

        srna = read_from_stack(stack=self.transform_stack, file='srna',
                               default_columns=self.srna_columns, reader=read_txt,
                               skiprows=38, names=self.srna_columns)

        sigma = read_from_stack(stack=self.transform_stack, file='sigma',
                                default_columns=self.sigma_columns, reader=read_txt,
                                skiprows=36, names=self.sigma_columns)

        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = self.select_columns(gene, 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                                   'ncbi_protein', 'genbank_accession', 'refseq_accession',
                                   'uniprot_accession', 'sequence', 'strand', 'start', 'stop',
                                   'gene_name_lower', 'gene_id', 'protrend_id')
        gene = gene.rename(columns={'protrend_id': 'gene_protrend_id'})

        tf = self._transform_tf(tf, gene)
        tf = tf.drop(columns=['site_length', 'symmetry', 'tf_internal_comment', 'key_id_org',
                              'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence'])
        srna = self._transform_srna(srna, gene)
        srna = srna.drop(columns=['srna_gene_regulated_id', 'srna_tu_regulated_id',
                                  'srna_function', 'srna_posleft', 'srna_posright', 'srna_sequence',
                                  'srna_regulation_type', 'srna_mechanis', 'srna_note'])
        sigma = self._transform_sigma(sigma, gene)
        sigma = sigma.drop(columns=['sigma_synonyms', 'sigma_coregulators', 'sigma_notes',
                                    'sigma_sigmulon_genes', 'key_id_org'])

        df = pd.concat([tf, srna, sigma], axis=0)
        self._stack_transformed_nodes(df)
        return df
