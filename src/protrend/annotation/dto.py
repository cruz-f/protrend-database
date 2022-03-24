from dataclasses import dataclass, fields, field
from typing import List, Any

from protrend.utils.set_list import set_list_field, SetList


@dataclass
class DTO:
    input_value: Any = field(default=None)

    def to_dict(self, ignore_fields: SetList[str] = None):

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
    locus_tag: SetList[str] = set_list_field(output='take_first')
    name: SetList[str] = set_list_field(output='take_first')
    synonyms: SetList[str] = set_list_field(output='take_all')
    function: SetList[str] = set_list_field(output='take_first')
    description: SetList[str] = set_list_field(output='take_first')
    ncbi_gene: SetList[str] = set_list_field(output='take_first')
    ncbi_protein: SetList[str] = set_list_field(output='take_first')
    genbank_accession: SetList[str] = set_list_field(output='take_first')
    refseq_accession: SetList[str] = set_list_field(output='take_first')
    uniprot_accession: SetList[str] = set_list_field(output='take_first')
    sequence: SetList[str] = set_list_field(output='take_first')
    strand: SetList[str] = set_list_field(output='take_first')
    start: List[int] = set_list_field(output='take_first')
    stop: List[int] = set_list_field(output='take_first')


@dataclass
class EffectorDTO(DTO):
    name: SetList[str] = set_list_field(output='take_first')
    synonyms: SetList[str] = set_list_field(output='take_all')
    mechanism: SetList[str] = set_list_field(output='take_first')
    kegg_compounds: SetList[str] = set_list_field(output='take_all')


@dataclass
class OrganismDTO(DTO):
    name: SetList[str] = set_list_field(output='take_first')
    species: SetList[str] = set_list_field(output='take_first')
    strain: SetList[str] = set_list_field(output='take_first')
    ncbi_taxonomy: SetList[str] = set_list_field(output='take_first')
    refseq_accession: SetList[str] = set_list_field(output='take_first')
    refseq_ftp: SetList[str] = set_list_field(output='take_first')
    genbank_accession: SetList[str] = set_list_field(output='take_first')
    genbank_ftp: SetList[str] = set_list_field(output='take_first')
    ncbi_assembly: SetList[str] = set_list_field(output='take_first')
    assembly_accession: SetList[str] = set_list_field(output='take_first')


@dataclass
class PathwayDTO(DTO):
    name: SetList[str] = set_list_field(output='take_first')
    synonyms: SetList[str] = set_list_field(output='take_all')
    kegg_pathways: SetList[str] = set_list_field(output='take_all')


@dataclass
class PublicationDTO(DTO):
    pmid: SetList[str] = set_list_field(output='take_first')
    doi: SetList[str] = set_list_field(output='take_first')
    title: SetList[str] = set_list_field(output='take_first')
    author: SetList[str] = set_list_field(output='take_first')
    year: SetList[str] = set_list_field(output='take_first')
