from typing import List, Union

import pandas as pd

from protrend.io import read_from_stack, read_json_lines, read_json_frame
from protrend.model import Gene
from protrend.annotation import annotate_genes, GeneDTO
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.collectf.regulator import RegulatorTransformer
from protrend.utils.processors import (take_first, flatten_set_list, to_int_str, apply_processors, rstrip, lstrip,
                                       to_list, to_set_list)
from protrend.utils import SetList


class GeneTransformer(CollectfTransformer,
                      source='collectf',
                      version='0.0.1',
                      node=Gene,
                      order=80,
                      register=True):
    default_transform_stack = {'gene': 'Gene.json', 'regulator': 'integrated_regulator.json'}
    columns = SetList(['name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession',
                       'uniprot_accession', 'sequence', 'strand', 'start', 'stop', 'regulon',
                       'operon', 'tfbs', 'regulator_protrend_id',
                       'regulator_uniprot_accession', 'ncbi_taxonomy', 'organism_protrend_id',
                       'locus_tag_old', 'locus_tag', 'protrend_id'])
    read_columns = SetList(['locus_tag', 'regulon', 'operon', 'tfbs'])

    def _transform_gene(self, gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        gene = apply_processors(gene, locus_tag=[rstrip, lstrip])
        gene = self.group_by(df=gene, column='locus_tag', aggregation={}, default=flatten_set_list)

        gene = apply_processors(gene, regulon=to_list)
        gene = gene.explode('regulon')

        gene = pd.merge(gene, regulator, left_on='regulon', right_on='regulator_uniprot_accession')
        gene = self.drop_duplicates(df=gene, subset=['locus_tag', 'organism_protrend_id'],
                                    perfect_match=True, preserve_nan=True)

        aggregation = {'regulon': to_set_list,
                       'regulator_uniprot_accession': to_set_list,
                       'regulator_protrend_id': to_set_list,
                       'ncbi_taxonomy': take_first,
                       'organism_protrend_id': take_first}
        gene = self.group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set_list)

        gene['locus_tag_old'] = gene['locus_tag']

        gene = self.create_input_value(df=gene, col='locus_tag')
        return gene

    @staticmethod
    def _annotate_genes(loci: List[Union[None, str]], taxa: List[str]):
        dtos = [GeneDTO(input_value=locus) for locus in loci]
        annotate_genes(dtos=dtos, loci=loci, taxa=taxa)

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

        genes = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (genes['strand'] != 'reverse') & (genes['strand'] != 'forward')
        genes.loc[strand_mask, 'strand'] = None
        return genes

    def transform(self):
        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=self.read_columns, reader=read_json_lines)
        non_empty_gene = gene['locus_tag'] != ''
        gene = gene[non_empty_gene]

        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'uniprot_accession',
                                        'ncbi_taxonomy', 'organism_protrend_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id',
                                              'uniprot_accession': 'regulator_uniprot_accession'})
        regulator = apply_processors(regulator, ncbi_taxonomy=to_int_str)

        gene = self._transform_gene(gene=gene, regulator=regulator)

        loci = gene['input_value'].tolist()
        taxa = gene['ncbi_taxonomy'].tolist()
        genes = self._annotate_genes(loci, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_collectf'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_collectf')

        df = df.drop(columns=['input_value'])

        df = apply_processors(df, ncbi_taxonomy=to_int_str)

        self._stack_transformed_nodes(df)

        return df
