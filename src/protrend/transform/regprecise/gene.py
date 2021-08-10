from collections import defaultdict
from typing import List, Union

import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.model.model import Source, Organism
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

    def _read_gene(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('gene')

        if file_path:
            df = read_json_lines(file_path)

        else:
            df = pd.DataFrame(columns=['locus_tag', 'name', 'function', 'url', 'regulon', 'operon', 'tfbs'])

        return df

    def _read_regulator(self) -> pd.DataFrame:
        file_path = self._transform_stack.get('regulator')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=RegulatorTransformer.columns)

        df = df.rename(columns={'protrend_id': 'regulator_protrend_id'})

        df = df[['organism_protrend_id', 'genome_id', 'ncbi_taxonomy',
                 'regulator_protrend_id', 'regulon_id']]

        return df

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
        gene = self._read_gene()
        regulator = self._read_regulator()

        gene = self._transform_gene(gene=gene, regulator=regulator)

        loci = gene['input_value'].tolist()
        names = gene['name'].tolist()
        taxa = gene['ncbi_taxonomy'].tolist()

        genes = self._annotate_genes(loci, names, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_regprecise'))

        # TODO: choose annotation if available

        df['locus_tag'] = df['locus_tag_annotation']
        df['name'] = df['name_annotation']
        df['function'] = df['function_annotation']

        df = df.drop(['input_value',
                      'name_annotation',
                      'name_regprecise',
                      'locus_tag_annotation',
                      'function_annotation',
                      'function_regprecise'], axis=1)

        if df.empty:
            df = self.make_empty_frame()

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df


class GeneToSourceConnector(DefaultConnector):
    default_settings = GeneToSource

    def _read_gene(self) -> pd.DataFrame:
        file_path = self._connect_stack.get('gene')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=GeneTransformer.columns)

        return df

    def _read_source(self) -> pd.DataFrame:
        file_path = self._connect_stack.get('source')

        if file_path:
            df = read_csv(file_path)

        else:
            df = pd.DataFrame(columns=SourceTransformer.columns)

        return df

    def connect(self):

        gene = self._read_gene()
        gene = gene.explode('regulon')
        source = self._read_source()

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
