from typing import Union

import pandas as pd

from protrend.annotation import (OrganismDTO, annotate_organisms)
from protrend.log import ProtrendLogger
from ._utils import get_values
from protrend.utils import SetList
from protrend.transform.transformer import Transformer


class OrganismMixIn:
    species = ['']
    strain = ['']
    ncbi_taxonomy = [None]
    refseq_accession = ['']
    refseq_ftp = ['']
    genbank_accession = ['']
    genbank_ftp = ['']
    ncbi_assembly = [None]
    assembly_accession = ['']
    name = ['']

    columns = SetList(['protrend_id', 'name', 'species', 'strain', 'ncbi_taxonomy', 'refseq_accession', 'refseq_ftp',
                       'genbank_accession', 'genbank_ftp', 'ncbi_assembly', 'assembly_accession'])

    @staticmethod
    def annotate_organisms(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        organisms = [OrganismDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(organisms)} organisms')

        identifiers = get_values(df, 'ncbi_taxonomy')
        names = get_values(df, 'name')

        iterator = zip(
            ('ncbi_taxonomy', 'name',),
            (identifiers, names)
        )

        params = [param for param, value in iterator if value is not None]
        params = ', '.join(params)

        ProtrendLogger.log.info(f'Annotating with the following params: {params}')

        annotate_organisms(dtos=organisms, identifiers=identifiers, names=names)

        organisms_dict = [dto.to_dict() for dto in organisms]

        organisms_df = pd.DataFrame(organisms_dict)
        return organisms_df

    def transform(self: Union[Transformer, 'OrganismMixIn']):
        org = dict(name=self.name,
                   species=self.species,
                   strain=self.strain,
                   ncbi_taxonomy=self.ncbi_taxonomy,
                   refseq_accession=self.refseq_accession,
                   refseq_ftp=self.refseq_ftp,
                   genbank_accession=self.genbank_accession,
                   genbank_ftp=self.genbank_ftp,
                   ncbi_assembly=self.ncbi_assembly,
                   assembly_accession=self.assembly_accession)

        df = pd.DataFrame(org, index=list(range(len(self.name))))

        self.stack_transformed_nodes(df)

        return df
