from typing import List

import pandas as pd

from protrend.bioapis import entrez_summary
from protrend.io import read_json_lines, read_from_stack
from protrend.model.model import Organism
from protrend.transform import OrganismDTO
from protrend.transform.annotation import annotate_organisms
from protrend.transform.collectf.base import CollectfTransformer
from protrend.transform.processors import apply_processors, rstrip, lstrip, to_int_str, take_last, flatten_set
from protrend.utils.miscellaneous import is_null


class OrganismTransformer(CollectfTransformer):
    default_node = Organism
    default_node_factors = ('ncbi_taxonomy', 'name')
    default_transform_stack = {'organism': 'Organism.json'}
    default_order = 100
    columns = {'protrend_id',
               'genome_accession', 'taxonomy', 'regulon', 'tfbs'
               'name', 'species', 'strain',
               'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
               'genbank_accession', 'genbank_ftp',
               'ncbi_assembly', 'assembly_accession'}
    read_columns = {'name', 'genome_accession', 'taxonomy', 'regulon', 'tfbs'}

    def _transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        aggregation = {'genome_accession': take_last, 'taxonomy': take_last}
        organism = self.group_by(df=organism, column='name', aggregation=aggregation, default=flatten_set)

        organism = apply_processors(organism, name=[rstrip, lstrip], genome_accession=[rstrip, lstrip],
                                    taxonomy=to_int_str)

        organism = self.create_input_value(organism, col='name')

        return organism

    @staticmethod
    def _transform_organisms(nucleotide: List[str], names: List[str]):

        identifiers = []
        for acc in nucleotide:

            ncbi_tax = None

            if not is_null(acc):
                summary = entrez_summary(identifier=acc, db='nucleotide')

                if 'TaxId' in summary:
                    ncbi_tax = str(summary['TaxId'])

            identifiers.append(ncbi_tax)

        dtos = [OrganismDTO(input_value=name) for name in names]
        annotate_organisms(dtos=dtos, identifiers=identifiers, names=names)

        # name: List[str]
        # species: List[str]
        # strain: List[str]
        # ncbi_taxonomy: List[int]
        # refseq_accession: List[str]
        # refseq_ftp: List[str]
        # genbank_accession: List[str]
        # genbank_ftp: List[str]
        # ncbi_assembly: List[str]
        # assembly_accession: List[str]

        return pd.DataFrame([dto.to_dict() for dto in dtos])

    def transform(self):
        organism = read_from_stack(stack=self.transform_stack, file='organism',
                                   default_columns=self.read_columns, reader=read_json_lines)
        organism = self._transform_organism(organism)

        names = organism['input_value'].tolist()
        nucleotide = organism['genome_accession'].tolist()
        organisms = self._transform_organisms(nucleotide, names)

        df = pd.merge(organisms, organism, on='input_value', suffixes=('_annotation', '_regprecise'))

        df = self.merge_columns(df=df, column='name', left='name_annotation', right='name_regprecise')

        df = df.drop(columns=['input_value'])

        df = apply_processors(df, ncbi_taxonomy=to_int_str, ncbi_assembly=to_int_str)

        self._stack_transformed_nodes(df)

        return df
