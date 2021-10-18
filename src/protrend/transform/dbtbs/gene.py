from typing import List

import pandas as pd

from protrend.io.json import read_json_lines
from protrend.io.utils import read_from_stack
from protrend.model.model import Gene
from protrend.transform.annotation import annotate_genes
from protrend.transform.dbtbs.base import DBTBSTransformer
from protrend.transform.dto import GeneDTO
from protrend.transform.processors import (rstrip, lstrip, apply_processors, lower_case)
from protrend.utils import SetList


class GeneTransformer(DBTBSTransformer):
    default_node = Gene
    default_transform_stack = {'gene': 'Gene.json'}
    default_order = 100
    columns = SetList(['protrend_id', 'locus_tag', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession',
                       'uniprot_accession', 'sequence', 'strand', 'start', 'stop',
                       'url', 'position', 'cog_id', 'conversed_groups', 'tf',
                       'operon', 'tfbs'])
    read_columns = SetList(['name', 'url', 'synonyms', 'strand', 'position', 'cog_id', 'conversed_groups', 'tf',
                            'operon', 'tfbs', 'function'])

    def _transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        gene = gene.drop(columns=['synonyms', 'strand', 'function'])

        gene = gene.explode(column='name')

        # filter nan and duplicates
        gene = self.drop_duplicates(df=gene, subset=['name'], perfect_match=True, preserve_nan=True)
        gene = gene.dropna(subset=['name'])

        gene = apply_processors(gene, name=[rstrip, lstrip])
        gene = self.create_input_value(df=gene, col='name')
        gene['input_value'] = apply_processors(df=gene, input_value=lower_case)

        return gene

    @staticmethod
    def _annotate_genes(names: List[str], taxa: List[str]):
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
        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=self.read_columns, reader=read_json_lines)

        gene = self._transform_gene(gene)

        names = gene['input_value'].tolist()
        taxa = gene['224308'] * len(names)

        genes = self._annotate_genes(names, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_regulondb'))

        df = df.dropna(subset=['locus_tag'])

        df['old_name'] = df['name_regulondb']
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regulondb')
        df = df.rename(columns={'old_name': 'name_regulondb'})

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)
        return df
