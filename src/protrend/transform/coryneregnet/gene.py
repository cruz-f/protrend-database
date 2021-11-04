from typing import List, Union

import pandas as pd

from protrend.model.model import Gene
from protrend.transform import GeneDTO
from protrend.annotation import annotate_genes
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.utils.processors import apply_processors, rstrip, lstrip
from protrend.utils import SetList


class GeneTransformer(CoryneRegNetTransformer):
    default_node = Gene
    default_transform_stack = {'bsub': 'bsub_regulation.csv',
                               'cglu': 'cglu_regulation.csv',
                               'ecol': 'ecol_regulation.csv',
                               'mtub': 'mtub_regulation.csv'}
    default_order = 100
    columns = SetList(['protrend_id', 'locus_tag', 'name', 'synonyms', 'function', 'description', 'ncbi_gene',
                       'ncbi_protein', 'genbank_accession', 'refseq_accession', 'uniprot_accession', 'sequence',
                       'strand', 'start', 'stop',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source', 'taxonomy'])

    def _transform_regulation(self, regulation: pd.DataFrame) -> pd.DataFrame:
        regulation = self.drop_duplicates(df=regulation, subset=['TG_locusTag'], perfect_match=True, preserve_nan=True)
        regulation = regulation.dropna(subset=['TG_locusTag'])

        regulation = apply_processors(regulation, TF_locusTag=[rstrip, lstrip], TF_name=[rstrip, lstrip])

        regulation['locus_tag'] = regulation['TG_locusTag']
        regulation['name'] = regulation['TG_name']

        regulation = self.create_input_value(df=regulation, col='locus_tag')
        return regulation

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
        regulation = self._build_regulations()

        regulation = self._transform_regulation(regulation)

        loci = regulation['input_value'].tolist()
        names = regulation['TF_name'].tolist()
        taxa = regulation['taxonomy'].tolist()

        genes = self._annotate_genes(loci, names, taxa)

        df = pd.merge(genes, regulation, on='input_value', suffixes=('_annotation', '_coryneregnet'))

        df = self.merge_columns(df=df, column='locus_tag', left='locus_tag_annotation', right='locus_tag_coryneregnet')
        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_coryneregnet')

        df = df.drop(columns=['input_value'])

        self._stack_transformed_nodes(df)

        return df
