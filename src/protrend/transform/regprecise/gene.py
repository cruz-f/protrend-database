from typing import List, Union

import pandas as pd

from protrend.io.json import read_json_lines, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.transform.annotation import annotate_genes
from protrend.transform.connector import DefaultConnector
from protrend.transform.transformer import Transformer
from protrend.transform.dto import GeneDTO
from protrend.transform.processors import rstrip, lstrip, apply_processors, take_last, flatten_set, to_list
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.regprecise.settings import GeneSettings, GeneToSource, GeneToOrganism


class GeneTransformer(Transformer):
    default_settings = GeneSettings
    columns = {'protrend_id',
               'locus_tag', 'name', 'synonyms', 'function', 'description',
               'ncbi_gene', 'ncbi_protein', 'genbank_accession',
               'refseq_accession', 'uniprot_accession',
               'sequence', 'strand', 'position_left', 'position_right',
               'organism_protrend_id', 'genome_id', 'ncbi_taxonomy',
               'regulator_protrend_id', 'regulon_id', 'locus_tag_old',
               'regulon', 'operon', 'tfbs'}

    read_columns = {'locus_tag', 'name', 'function', 'url', 'regulon', 'operon', 'tfbs'}

    def _transform_gene(self, gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        apply_processors(rstrip, lstrip, df=gene, col='locus_tag')
        apply_processors(rstrip, lstrip, df=gene, col='name')

        aggregation = {'name': take_last, 'function': take_last, 'url': set}
        gene = self.group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set)

        apply_processors(to_list, df=gene, col='regulon')
        gene = gene.explode('regulon')

        gene = pd.merge(gene, regulator, how='left', left_on='regulon', right_on='regulon_id')

        aggregation = {'name': take_last, 'function': take_last,
                       'organism_protrend_id': take_last, 'genome_id': take_last, 'ncbi_taxonomy': take_last,
                       'regulator_protrend_id': take_last, 'regulon_id': take_last, 'regulon': set}
        gene = self.group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set)

        gene = self.create_input_value(df=gene, col='locus_tag')
        return gene

    @staticmethod
    def _annotate_genes(loci: List[Union[None, str]], names: List[str], taxa: List[int]):
        dtos = [GeneDTO(input_value=locus) for locus in loci]
        annotate_genes(dtos=dtos, loci=loci, names=names, taxa=taxa)

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

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        gene = read_from_stack(stack=self._transform_stack, file='gene',
                               default_columns=self.read_columns, reader=read_json_lines)

        regulator = read_from_stack(stack=self._transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'genome_id', 'ncbi_taxonomy', 'regulon_id',
                                        'organism_protrend_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})

        gene = self._transform_gene(gene=gene, regulator=regulator)
        gene['locus_tag_old'] = gene['locus_tag']

        loci = gene['input_value'].tolist()
        names = gene['name'].tolist()
        taxa = gene['ncbi_taxonomy'].tolist()
        genes = self._annotate_genes(loci, names, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_regprecise')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')
        df = self.merge_columns(df=df, column='function', left='function_annotation', right='function_regprecise')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df


class GeneToSourceConnector(DefaultConnector):
    default_settings = GeneToSource

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = gene.explode('regulon')
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        from_identifiers = gene['protrend_id'].tolist()
        size = len(from_identifiers)

        protrend_id = source['protrend_id'].iloc[0]
        to_identifiers = [protrend_id] * size

        kwargs = dict(url=gene['url'].tolist(),
                      external_identifier=gene['regulon'].tolist(),
                      key=['regulon_id'] * size)

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers,
                                  kwargs=kwargs)

        self.stack_csv(df)


class GeneToOrganismConnector(DefaultConnector):
    default_settings = GeneToOrganism

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)

        from_identifiers = gene['protrend_id'].tolist()
        to_identifiers = gene['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)
