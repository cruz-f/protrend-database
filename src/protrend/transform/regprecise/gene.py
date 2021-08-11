from typing import List, Union

import pandas as pd

from protrend.io.utils import read_from_stack
from protrend.transform.annotation.gene import annotate_genes
from protrend.transform.connector import DefaultConnector
from protrend.transform.dto import GeneDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, take_first
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.settings import GeneSettings, GeneToSource
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import DefaultTransformer
from protrend.utils.miscellaneous import take_last, flatten_list


class GeneTransformer(DefaultTransformer):
    default_settings = GeneSettings
    columns = {'protrend_id',
               'locus_tag', 'name', 'synonyms', 'function', 'description',
               'ncbi_gene', 'ncbi_protein', 'genbank_accession',
               'refseq_accession', 'uniprot_accession',
               'sequence', 'strand', 'position_left', 'position_right',
               'annotation_score',
               'organism_protrend_id', 'genome_id', 'ncbi_taxonomy',
               'regulator_protrend_id', 'regulon_id', 'locus_tag_regprecise'}
    read_columns = {'locus_tag', 'name', 'function', 'url', 'regulon', 'operon', 'tfbs'}

    @staticmethod
    def _transform_gene(gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:

        apply_processors(rstrip,
                         lstrip,
                         df=gene,
                         col='locus_tag')

        apply_processors(rstrip,
                         lstrip,
                         df=gene,
                         col='name')

        aggregation_functions = {'name': take_last,
                                 'function': take_last,
                                 'url': set,
                                 'regulon': flatten_list,
                                 'operon': flatten_list,
                                 'tfbs': flatten_list}

        gene = gene.groupby(gene['locus_tag']).aggregate(aggregation_functions)
        gene = gene.reset_index()

        gene['regulon_id'] = gene['regulon']

        # keeping only one, since we only want to get the ncbi taxonomy of each gene.
        apply_processors(take_first,
                         df=gene,
                         col='regulon_id')

        df = pd.merge(gene, regulator, on='regulon_id')

        df['input_value'] = df['locus_tag']

        return df

    @staticmethod
    def _annotate_genes(loci: List[Union[None, str]], names: List[str], taxa: List[int]):

        dtos = [GeneDTO(input_value=locus) for locus in loci]
        annotate_genes(dtos=dtos, loci=loci, names=names, taxa=taxa)

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
        # position_left: List[int]
        # position_right: List[int]
        # annotation_score: int

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        gene = read_from_stack(tl=self, file='gene', json=True, default_columns=self.read_columns)

        regulator = read_from_stack(tl=self, file='regulator', json=False, default_columns=RegulatorTransformer.columns)
        regulator = regulator[['protrend_id', 'genome_id', 'ncbi_taxonomy', 'regulator_protrend_id', 'regulon_id']]
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        gene = self._transform_gene(gene=gene, regulator=regulator)

        loci = gene['input_value'].tolist()
        names = gene['name'].tolist()
        taxa = gene['ncbi_taxonomy'].tolist()

        genes = self._annotate_genes(loci, names, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_regprecise'))

        old_loci = df['locus_tag_regprecise'].tolist()

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation',
                                right='locus_tag_regprecise', fill='')
        df = self.merge_columns(df=df, column='name', left='name_annotation',
                                right='name_regprecise', fill='')
        df = self.merge_columns(df=df, column='function', left='function_annotation',
                                right='function_regprecise', fill='')

        df['locus_tag_regprecise'] = old_loci

        df = df.drop(['input_value'], axis=1)

        if df.empty:
            df = self.make_empty_frame()

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df


class GeneToSourceConnector(DefaultConnector):
    default_settings = GeneToSource

    def connect(self):

        gene = read_from_stack(tl=self, file='gene', json=False, default_columns=GeneTransformer.columns)
        gene = gene.explode('regulon')
        source = read_from_stack(tl=self, file='source', json=False, default_columns=SourceTransformer.columns)

        from_identifiers = gene['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=gene['url'].tolist(),
                      external_identifier=gene['regulon'].tolist(),
                      key=['regulon_id'] * size)

        df = self.make_connection(size=size,
                                  from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)
