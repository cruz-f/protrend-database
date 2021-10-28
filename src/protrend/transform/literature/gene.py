from typing import List, Union, Tuple

import pandas as pd

from protrend.model.model import Gene
from protrend.transform import GeneDTO
from protrend.transform.annotation import annotate_genes
from protrend.transform.literature.base import LiteratureTransformer
from protrend.transform.processors import apply_processors, rstrip, lstrip, to_set_list
from protrend.utils import SetList


class GeneTransformer(LiteratureTransformer):
    default_node = Gene
    default_order = 100
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession', 'sequence',
                       'strand', 'start', 'stop',
                       'regulator_locus_tag', 'regulator_name', 'operon', 'genes_locus_tag',
                       'genes_name', 'regulatory_effect', 'evidence', 'effector', 'mechanism',
                       'publication', 'taxonomy', 'source', 'network_id'])

    @staticmethod
    def _explode_loci_names(network):
        network = apply_processors(network, genes_locus_tag=to_set_list, genes_name=to_set_list)

        def reshape_genes_locus_tag(genes_locus_tag: list, genes_name: list) -> Tuple[list, list]:

            if len(genes_locus_tag) != len(genes_name):

                if not genes_locus_tag:
                    return genes_name

            return genes_locus_tag

        def reshape_genes_names(genes_locus_tag: list, genes_name: list) -> Tuple[list, list]:

            if len(genes_locus_tag) != len(genes_name):

                if not genes_locus_tag:
                    return genes_name

            return genes_name[:len(genes_locus_tag)]

        network = apply_processors(network, genes_locus_tag=reshape_genes_locus_tag, genes_name=reshape_genes_names)

        network = network.set_index(['operon']).apply(pd.Series.explode).reset_index()
        return network

    def _transform_gene(self, network: pd.DataFrame) -> pd.DataFrame:
        network = self._explode_loci_names(network)

        network = apply_processors(network, genes_locus_tag=[rstrip, lstrip], genes_name=[rstrip, lstrip])

        network = self.drop_duplicates(df=network, subset=['genes_locus_tag'], perfect_match=True, preserve_nan=True)
        network = network.dropna(subset=['genes_locus_tag'])

        network['locus_tag'] = network['genes_locus_tag']
        network['name'] = network['genes_name']

        network = self.create_input_value(df=network, col='locus_tag')
        return network

    @staticmethod
    def _annotate_genes(loci: List[Union[None, str]], names: List[str], taxa: List[str]):
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
        # start: List[int]
        # stop: List[int]

        genes = pd.DataFrame([dto.to_dict() for dto in dtos])
        strand_mask = (genes['strand'] != 'reverse') & (genes['strand'] != 'forward')
        genes.loc[strand_mask, 'strand'] = None
        return genes

    def transform(self):
        network = self._build_network()
        gene = self._transform_gene(network)

        loci = gene['input_value'].tolist()
        names = gene['genes_name'].tolist()
        taxa = gene['taxonomy'].tolist()

        genes = self._annotate_genes(loci, names, taxa)

        df = pd.merge(genes, gene, on='input_value', suffixes=('_annotation', '_literature'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_literature')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_literature')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df
