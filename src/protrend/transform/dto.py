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
    locus_tag: List[str] = set_list_field(output='take_first', init=False)
    name: List[str] = set_list_field(output='take_first', init=False)
    synonyms: List[str] = set_list_field(output='take_all', init=False)
    function: List[str] = set_list_field(output='take_first', init=False)
    description: List[str] = set_list_field(output='take_first', init=False)
    ncbi_gene: List[str] = set_list_field(output='take_first', init=False)
    ncbi_protein: List[str] = set_list_field(output='take_first', init=False)
    genbank_accession: List[str] = set_list_field(output='take_first', init=False)
    refseq_accession: List[str] = set_list_field(output='take_first', init=False)
    uniprot_accession: List[str] = set_list_field(output='take_first', init=False)
    sequence: List[str] = set_list_field(output='take_first', init=False)
    strand: List[str] = set_list_field(output='take_first', init=False)
    start: List[int] = set_list_field(output='take_first', init=False)
    stop: List[int] = set_list_field(output='take_first', init=False)


@dataclass
class EffectorDTO(DTO):
    name: List[str] = set_list_field(output='take_first', init=False)
    synonyms: List[str] = set_list_field(output='take_all', init=False)
    mechanism: List[str] = set_list_field(output='take_first', init=False)
    kegg_compounds: List[str] = set_list_field(output='take_all', init=False)


@dataclass
class OrganismDTO(DTO):
    name: List[str] = set_list_field(output='take_first', init=False)
    species: List[str] = set_list_field(output='take_first', init=False)
    strain: List[str] = set_list_field(output='take_first', init=False)
    family: List[str] = set_list_field(output='take_first', init=False)
    phylum: List[str] = set_list_field(output='take_first', init=False)
    ncbi_taxonomy: List[str] = set_list_field(output='take_first', init=False)
    refseq_accession: List[str] = set_list_field(output='take_first', init=False)
    refseq_ftp: List[str] = set_list_field(output='take_first', init=False)
    genbank_accession: List[str] = set_list_field(output='take_first', init=False)
    genbank_ftp: List[str] = set_list_field(output='take_first', init=False)
    ncbi_assembly: List[str] = set_list_field(output='take_first', init=False)
    assembly_accession: List[str] = set_list_field(output='take_first', init=False)


@dataclass
class PathwayDTO(DTO):
    name: List[str] = set_list_field(output='take_first', init=False)
    synonyms: List[str] = set_list_field(output='take_all', init=False)
    kegg_pathways: List[str] = set_list_field(output='take_all', init=False)


@dataclass
class PublicationDTO(DTO):
    pmid: List[str] = set_list_field(output='take_first', init=False)
    doi: List[str] = set_list_field(output='take_first', init=False)
    title: List[str] = set_list_field(output='take_first', init=False)
    author: List[str] = set_list_field(output='take_first', init=False)
    year: List[str] = set_list_field(output='take_first', init=False)
