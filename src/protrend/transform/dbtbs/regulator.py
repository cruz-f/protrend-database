import pandas as pd
from Bio import SeqIO

from protrend.io import read_json_lines, read_from_stack
from protrend.model import Regulator
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.mix_ins import SequenceMixIn, GeneMixIn
from protrend.transform.transformations import drop_empty_string, drop_duplicates, create_input_value
from protrend.utils import SetList
from protrend.utils.processors import rstrip, lstrip, apply_processors, to_list_nan, take_first


class RegulatorTransformer(GeneMixIn, SequenceMixIn, DBTBSTransformer,
                           source='dbtbs',
                           version='0.0.4',
                           node=Regulator,
                           order=100,
                           register=True):
    default_transform_stack = {'tf': 'TranscriptionFactor.json', 'sequence': 'sequence.gb'}
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'mechanism',
                       'family', 'domain', 'domain_description', 'description', 'url',
                       'type', 'consensus_sequence', 'comment', 'subti_list', 'gene', 'tfbs',
                       'name_lower', 'name_dbts'])

    read_columns = SetList(['name', 'family', 'domain', 'domain_description', 'description', 'url',
                            'type', 'consensus_sequence', 'comment', 'subti_list', 'gene', 'tfbs'])

    @staticmethod
    def transform_tf(tf: pd.DataFrame, sequence: pd.DataFrame) -> pd.DataFrame:
        tf = tf.assign(name_dbtbs=tf['name'].copy())

        tf = apply_processors(tf, name=[rstrip, lstrip], family=[to_list_nan, take_first, rstrip, lstrip])

        # filter nan and duplicates
        tf = tf.dropna(subset=['name'])
        tf = drop_empty_string(tf, 'name')
        tf = drop_duplicates(df=tf, subset=['name'])

        tf = tf.assign(name_lower=tf['name'].str.lower(), mechanism=None)

        # sigma factors
        tf_sigma_mask = tf['family'] == 'Sigma factors'
        tf.loc[tf_sigma_mask, 'mechanism'] = 'sigma factor'

        # tf mechanism
        tf_not_sigma_mask = tf['family'] != 'Sigma factors'
        tf.loc[tf_not_sigma_mask, 'mechanism'] = 'transcription factor'

        tf = pd.merge(tf, sequence, on='name_lower')

        tf = create_input_value(df=tf, col='locus_tag')
        return tf

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, key='tf',
                             columns=self.read_columns, reader=read_json_lines)

        gb_file = self.transform_stack['sequence']
        sequence = SeqIO.read(gb_file, "genbank")
        sequence = self.transform_sequence(sequence)

        regulators = self.transform_tf(tf=tf, sequence=sequence)
        annotated_regulators = self.annotate_genes(regulators)

        df = self.merge_annotations(annotated_regulators, regulators)

        df = df.drop(columns=['input_value'])

        self.stack_transformed_nodes(df)
        return df
