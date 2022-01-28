from dataclasses import dataclass, fields, field
from typing import List, Any

from protrend.utils.set_list import set_list_field


@dataclass
class DTO:
    input_value: Any = field(default=None)

    def to_dict(self, ignore_fields: List[str] = None):

        if not ignore_fields:
            ignore_fields = []

        result = []

        for f in fields(self):

            if f.name in ignore_fields:
                continue

            elif f.name == 'input_value':
                value = getattr(self, f.name)

            else:

                # ouput processors of the SetList are only used when the set list is called as function
                set_list = getattr(self, f.name)
                value = set_list()

            result.append((f.name, value))

        return dict(result)


@dataclass
class GeneDTO(DTO):
    locus_tag: List[str] = set_list_field(output='take_first')
    name: List[str] = set_list_field(output='take_first')
    synonyms: List[str] = set_list_field(output='take_all')
    function: List[str] = set_list_field(output='take_first')
    description: List[str] = set_list_field(output='take_first')
    ncbi_gene: List[str] = set_list_field(output='take_first')
    ncbi_protein: List[str] = set_list_field(output='take_first')
    genbank_accession: List[str] = set_list_field(output='take_first')
    refseq_accession: List[str] = set_list_field(output='take_first')
    uniprot_accession: List[str] = set_list_field(output='take_first')
    sequence: List[str] = set_list_field(output='take_first')
    strand: List[str] = set_list_field(output='take_first')
    start: List[int] = set_list_field(output='take_first')
    stop: List[int] = set_list_field(output='take_first')


@dataclass
class EffectorDTO(DTO):
    name: List[str] = set_list_field(output='take_first')
    synonyms: List[str] = set_list_field(output='take_all')
    mechanism: List[str] = set_list_field(output='take_first')
    kegg_compounds: List[str] = set_list_field(output='take_all')


@dataclass
class OrganismDTO(DTO):
    name: List[str] = set_list_field(output='take_first')
    species: List[str] = set_list_field(output='take_first')
    strain: List[str] = set_list_field(output='take_first')
    ncbi_taxonomy: List[str] = set_list_field(output='take_first')
    refseq_accession: List[str] = set_list_field(output='take_first')
    refseq_ftp: List[str] = set_list_field(output='take_first')
    genbank_accession: List[str] = set_list_field(output='take_first')
    genbank_ftp: List[str] = set_list_field(output='take_first')
    ncbi_assembly: List[str] = set_list_field(output='take_first')
    assembly_accession: List[str] = set_list_field(output='take_first')


@dataclass
class PathwayDTO(DTO):
    name: List[str] = set_list_field(output='take_first')
    synonyms: List[str] = set_list_field(output='take_all')
    kegg_pathways: List[str] = set_list_field(output='take_all')


@dataclass
class PublicationDTO(DTO):
    pmid: List[str] = set_list_field(output='take_first')
    doi: List[str] = set_list_field(output='take_first')
    title: List[str] = set_list_field(output='take_first')
    author: List[str] = set_list_field(output='take_first')
    year: List[str] = set_list_field(output='take_first')
