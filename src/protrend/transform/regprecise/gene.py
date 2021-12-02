from typing import List, Union

import pandas as pd

from protrend.io import read_json_lines, read_json_frame, read_from_stack
from protrend.model import Gene, Source, Organism
from protrend.annotation import annotate_genes, GeneDTO
from protrend.utils.processors import (rstrip, lstrip, apply_processors, take_last,
                                       flatten_set_list, to_list, to_int_str)
from protrend.transform.regprecise.base import RegPreciseTransformer, RegPreciseConnector
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.source import SourceTransformer
from protrend.utils import SetList


class GeneTransformer(RegPreciseTransformer,
                      source='regprecise',
                      version='0.0.0',
                      node=Gene,
                      order=80,
                      register=True):
    default_transform_stack = {'gene': 'Gene.json', 'regulator': 'integrated_regulator.json'}
    columns = SetList(['synonyms', 'description', 'ncbi_gene', 'ncbi_protein',
                       'genbank_accession', 'refseq_accession', 'uniprot_accession',
                       'sequence', 'strand', 'start', 'stop', 'url', 'regulon', 'operon',
                       'tfbs', 'regulator_protrend_id', 'genome_id', 'ncbi_taxonomy',
                       'regulon_id', 'organism_protrend_id', 'locus_tag_old', 'locus_tag',
                       'name', 'function', 'protrend_id'])
    read_columns = SetList(['locus_tag', 'name', 'function', 'url', 'regulon', 'operon', 'tfbs'])

    def _transform_gene(self, gene: pd.DataFrame, regulator: pd.DataFrame) -> pd.DataFrame:
        gene = apply_processors(gene, locus_tag=[rstrip, lstrip], name=[rstrip, lstrip])

        aggregation = {'name': take_last, 'function': take_last, 'url': set}
        gene = self.group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set_list)

        gene = apply_processors(gene, regulon=to_list)
        gene = gene.explode('regulon')
        gene = apply_processors(gene, regulon=to_int_str)

        gene = pd.merge(gene, regulator, how='left', left_on='regulon', right_on='regulon_id')

        aggregation = {'name': take_last, 'function': take_last,
                       'organism_protrend_id': take_last, 'genome_id': take_last, 'ncbi_taxonomy': take_last,
                       'regulator_protrend_id': take_last, 'regulon_id': take_last, 'regulon': set}
        gene = self.group_by(df=gene, column='locus_tag', aggregation=aggregation, default=flatten_set_list)

        gene['locus_tag_old'] = gene['locus_tag']

        gene = self.create_input_value(df=gene, col='locus_tag')
        return gene

    @staticmethod
    def _annotate_genes(loci: List[Union[None, str]], names: List[str], taxa: List[str]):
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
        # start: List[int]
        # stop: List[int]

        genes = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (genes['strand'] != 'reverse') & (genes['strand'] != 'forward')
        genes.loc[strand_mask, 'strand'] = None
        return genes

    def transform(self):
        gene = read_from_stack(stack=self.transform_stack, file='gene',
                               default_columns=self.read_columns, reader=read_json_lines)

        regulator = read_from_stack(stack=self.transform_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)
        regulator = self.select_columns(regulator, 'protrend_id', 'genome_id', 'ncbi_taxonomy', 'regulon_id',
                                        'organism_protrend_id')
        regulator = regulator.rename(columns={'protrend_id': 'regulator_protrend_id'})
        regulator = apply_processors(regulator, genome_id=to_int_str, ncbi_taxonomy=to_int_str, regulon_id=to_int_str)

        gene = self._transform_gene(gene=gene, regulator=regulator)

        loci = gene['input_value'].tolist()
        names = gene['name'].tolist()
        taxa = gene['ncbi_taxonomy'].tolist()
        genes = self._annotate_genes(loci, names, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_regprecise')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')
        df = self.merge_columns(df=df, column='function', left='function_annotation', right='function_regprecise')

        df = df.drop(columns=['input_value'])

        df = apply_processors(df, genome_id=to_int_str, ncbi_taxonomy=to_int_str, regulon_id=to_int_str)

        self.stack_transformed_nodes(df)

        return df


class GeneToSourceConnector(RegPreciseConnector,
                            source='regprecise',
                            version='0.0.0',
                            from_node=Gene,
                            to_node=Source,
                            register=True):
    default_connect_stack = {'gene': 'integrated_gene.json', 'source': 'integrated_source.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)
        gene = apply_processors(gene, regulon=to_list)
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

        self.stack_json(df)


class GeneToOrganismConnector(RegPreciseConnector,
                              source='regprecise',
                              version='0.0.0',
                              from_node=Gene,
                              to_node=Organism,
                              register=True):
    default_connect_stack = {'gene': 'integrated_gene.json'}

    def connect(self):
        gene = read_from_stack(stack=self._connect_stack, file='gene',
                               default_columns=GeneTransformer.columns, reader=read_json_frame)

        from_identifiers = gene['protrend_id'].tolist()
        to_identifiers = gene['organism_protrend_id'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_json(df)
