from typing import List

import pandas as pd

from protrend.io.json import read_json_lines
from protrend.io.utils import read_from_stack
from protrend.model.model import Regulator
from protrend.transform.annotation import annotate_genes
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.dto import GeneDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, lower_case
from protrend.utils import SetList


class RegulatorTransformer(DBTBSTransformer):
    default_node = Regulator
    default_transform_stack = {'tf': 'TranscriptionFactor.json'}
    default_order = 100
    columns = SetList(['locus_tag', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession',
                       'uniprot_accession', 'sequence', 'strand', 'start', 'stop', 'mechanism', 'protrend_id',
                       'name_regulondb', 'family', 'domain', 'domain_description', 'url',
                       'type', 'comment', 'operon', 'subti_list', 'consensus_sequence'])

    read_columns = SetList(['name', 'family', 'domain', 'domain_description', 'description', 'url',
                            'type', 'comment', 'operon', 'subti_list', 'consensus_sequence'])

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        tf = tf.drop(columns=['description'])

        # filter nan and duplicates
        tf = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True,
                                  preserve_nan=True)
        tf = tf.dropna(subset=['name'])

        tf['mechanism'] = 'transcription factor'

        tf = apply_processors(tf, name=[rstrip, lstrip])
        tf = self.create_input_value(df=tf, col='name')
        tf['input_value'] = apply_processors(df=tf, input_value=lower_case)

        return tf

    @staticmethod
    def _annotate_tfs(names: List[str], taxa: List[str]):
        dtos = [GeneDTO(input_value=name) for name in names]
        annotate_genes(dtos=dtos, names=names, taxa=taxa)

        for dto, name in zip(dtos, names):
            dto.synonyms.append(name)

        # locus_tag: List[str]
        # name: List[str]
        # synonyms: List[str]
        # function: List[str]
        # description: List[str]
        # ncbi_gene: List[str]
        # ncbi_protein: List[str]
        # genbank_accession: List[str]
        # refseq_accession: List[str]
        # uniprot_accession: List[str]
        # sequence: List[str]
        # strand: List[str]
        # start: List[int]
        # stop: List[int]

        df = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (df['strand'] != 'reverse') & (df['strand'] != 'forward')
        df.loc[strand_mask, 'strand'] = None
        return df

    def transform(self):
        tf = read_from_stack(stack=self.transform_stack, file='tf',
                             default_columns=self.read_columns, reader=read_json_lines)

        tf = self._transform_tf(tf)

        names = tf['input_value'].tolist()
        taxa = tf['224308'] * len(names)

        tfs = self._annotate_tfs(names, taxa)

        df = pd.merge(tfs, tf, on='input_value', suffixes=('_annotation', '_regulondb'))

        df = df.dropna(subset=['locus_tag'])

        df['old_name'] = df['name_regulondb']
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regulondb')
        df = df.rename(columns={'old_name': 'name_regulondb'})

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)
        return df
