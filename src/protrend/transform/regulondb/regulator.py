import pandas as pd

from protrend.io import read_from_stack, read_json_frame
from protrend.model import Regulator
from protrend.transform.regulondb.base import RegulondbTransformer, regulondb_reader
from protrend.transform.regulondb.gene import GeneTransformer
from protrend.transform.transformations import select_columns, drop_empty_string, drop_duplicates
from protrend.utils import SetList
from protrend.utils.processors import apply_processors, rstrip, lstrip


class RegulatorTransformer(RegulondbTransformer,
                           source='regulondb',
                           version='0.0.0',
                           node=Regulator,
                           order=90,
                           register=True):
    default_transform_stack = {'tf': 'transcription_factor.txt',
                               'srna': 'srna_interaction.txt',
                               'sigma': 'sigma_tmp.txt',
                               'gene': 'integrated_gene.json'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
                       'name_lower', 'gene_id', 'regulator_id',
                       'transcription_factor_id', 'transcription_factor_name', 'transcription_factor_family',
                       'srna_id', 'srna_gene_id',
                       'sigma_id', 'sigma_gene_id'])

    tf_columns = SetList(['transcription_factor_id', 'transcription_factor_name', 'site_length', 'symmetry',
                          'transcription_factor_family', 'tf_internal_comment', 'key_id_org',
                          'transcription_factor_note', 'connectivity_class', 'sensing_class', 'consensus_sequence'])
    srna_columns = SetList(['srna_id', 'srna_gene_id', 'srna_gene_regulated_id', 'srna_tu_regulated_id',
                            'srna_function', 'srna_posleft', 'srna_posright', 'srna_sequence',
                            'srna_regulation_type', 'srna_mechanis', 'srna_note'])
    sigma_columns = SetList(['sigma_id', 'sigma_name', 'sigma_synonyms', 'sigma_gene_id', 'sigma_gene_name',
                             'sigma_coregulators', 'sigma_notes', 'sigma_sigmulon_genes', 'key_id_org'])

    @staticmethod
    def transform_tf(tf: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        tf = select_columns(tf,
                            'transcription_factor_id', 'transcription_factor_name',
                            'transcription_factor_family')

        tf = tf.assign(name_lower=tf['transcription_factor_name'].str.lower(),
                       mechanism='transcription factor',
                       regulator_id=tf['transcription_factor_id'].copy())

        tf = apply_processors(tf, name_lower=[rstrip, lstrip])
        tf = tf.dropna(subset=['name_lower'])
        tf = drop_empty_string(tf, 'name_lower')
        tf = drop_duplicates(df=tf, subset=['name_lower'])

        tf = pd.merge(tf, gene, on='name_lower')
        return tf

    @staticmethod
    def transform_sigma(sigma: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        sigma = select_columns(sigma, 'sigma_id', 'sigma_gene_id')

        sigma = sigma.assign(regulator_id=sigma['sigma_id'].copy(),
                             mechanism='sigma factor')

        sigma = apply_processors(sigma, sigma_gene_id=[rstrip, lstrip])
        sigma = sigma.dropna(subset=['sigma_gene_id'])
        sigma = drop_empty_string(sigma, 'sigma_gene_id')
        sigma = drop_duplicates(df=sigma, subset=['sigma_gene_id'])

        sigma = pd.merge(sigma, gene, left_on='sigma_gene_id', right_on='gene_id')
        return sigma

    @staticmethod
    def transform_srna(srna: pd.DataFrame, gene: pd.DataFrame) -> pd.DataFrame:
        srna = select_columns(srna, 'srna_id', 'srna_gene_id')

        srna = srna.assign(regulator_id=srna['srna_gene_id'].copy(),
                           mechanism='small RNA (sRNA)')

        srna = apply_processors(srna, srna_gene_id=[rstrip, lstrip])
        srna = srna.dropna(subset=['srna_gene_id'])
        srna = drop_empty_string(srna, 'srna_gene_id')
        srna = drop_duplicates(df=srna, subset=['srna_gene_id'])

        srna = pd.merge(srna, gene, left_on='srna_gene_id', right_on='gene_id')
        return srna

    def transform(self):
        tf_reader = regulondb_reader(skiprows=38, names=self.tf_columns)
        tf = read_from_stack(stack=self.transform_stack, key='tf',
                             columns=self.tf_columns, reader=tf_reader)

        srna_reader = regulondb_reader(skiprows=38, names=self.srna_columns)
        srna = read_from_stack(stack=self.transform_stack, key='srna',
                               columns=self.srna_columns, reader=srna_reader)

        sigma_reader = regulondb_reader(skiprows=36, names=self.sigma_columns)
        sigma = read_from_stack(stack=self.transform_stack, key='sigma',
                                columns=self.sigma_columns, reader=sigma_reader)

        gene = read_from_stack(stack=self.transform_stack, key='gene',
                               columns=GeneTransformer.columns, reader=read_json_frame)

        gene = select_columns(gene, 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                              'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                              'sequence', 'strand', 'start', 'stop',
                              'name_lower', 'gene_id')

        tf = self.transform_tf(tf, gene)
        srna = self.transform_srna(srna, gene)
        sigma = self.transform_sigma(sigma, gene)

        df = pd.concat([tf, srna, sigma])
        df = df.reset_index(drop=True)
        self.stack_transformed_nodes(df)
        return df
