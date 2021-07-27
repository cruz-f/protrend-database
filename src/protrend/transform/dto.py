from dataclasses import dataclass, fields
from typing import List

from protrend.utils.unique_list import unique_field


@dataclass
class DTO:

    def to_dict(self, ignore_fields: List[str] = None):

        if not ignore_fields:
            ignore_fields = []

        result = []

        for f in fields(self):

            if f.name in ignore_fields:
                continue

            unique_list = getattr(self, f.name)
            value = unique_list()

            result.append((f.name, value))

        return dict(result)


@dataclass
class GeneDTO(DTO):
    locus_tag: List[str] = unique_field(output='take_first', init=False)
    name: List[str] = unique_field(output='take_first', init=False)
    synonyms: List[str] = unique_field(output='take_all', init=False)
    function: List[str] = unique_field(output='take_first', init=False)
    description: List[str] = unique_field(output='take_first', init=False)
    ncbi_gene: List[str] = unique_field(output='take_first', init=False)
    ncbi_protein: List[str] = unique_field(output='take_first', init=False)
    genbank_accession: List[str] = unique_field(output='take_first', init=False)
    refseq_accession: List[str] = unique_field(output='take_first', init=False)
    uniprot_accession: List[str] = unique_field(output='take_first', init=False)
    sequence: List[str] = unique_field(output='take_first', init=False)
    strand: List[str] = unique_field(output='take_first', init=False)
    position_left: List[int] = unique_field(output='take_first', init=False)
    position_right: List[int] = unique_field(output='take_first', init=False)
    annotation_score: int = 0


@dataclass
class EffectorDTO(DTO):
    name: List[str] = unique_field(output='take_first', init=False)
    synonyms: List[str] = unique_field(output='take_all', init=False)
    mechanism: List[str] = unique_field(output='take_first', init=False)
    kegg_compounds: List[str] = unique_field(output='take_all', init=False)


@dataclass
class OrganismDTO(DTO):
    name: List[str] = unique_field(output='take_first', init=False)
    species: List[str] = unique_field(output='take_first', init=False)
    strain: List[str] = unique_field(output='take_first', init=False)
    family: List[str] = unique_field(output='take_first', init=False)
    phylum: List[str] = unique_field(output='take_first', init=False)
    ncbi_taxonomy: List[str] = unique_field(output='take_first', init=False)
    refseq_accession: List[str] = unique_field(output='take_first', init=False)
    refseq_ftp: List[str] = unique_field(output='take_first', init=False)
    genbank_accession: List[str] = unique_field(output='take_first', init=False)
    genbank_ftp: List[str] = unique_field(output='take_first', init=False)
    ncbi_assembly: List[str] = unique_field(output='take_first', init=False)
    assembly_accession: List[str] = unique_field(output='take_first', init=False)


@dataclass
class PathwayDTO(DTO):
    name: List[str] = unique_field(output='take_first', init=False)
    synonyms: List[str] = unique_field(output='take_all', init=False)
    kegg_pathways: List[str] = unique_field(output='take_all', init=False)


@dataclass
class PublicationDTO(DTO):
    pmid: List[str] = unique_field(output='take_first', init=False)
    doi: List[str] = unique_field(output='take_first', init=False)
    title: List[str] = unique_field(output='take_first', init=False)
    author: List[str] = unique_field(output='take_first', init=False)
    year: List[str] = unique_field(output='take_first', init=False)
