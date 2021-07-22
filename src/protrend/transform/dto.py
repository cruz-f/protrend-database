from dataclasses import dataclass, field
from typing import Set


@dataclass
class GeneDTO:
    locus_tag: Set[str] = field(default_factory=set)
    name: Set[str] = field(default_factory=set)
    synonyms: Set[str] = field(default_factory=set)
    function: Set[str] = field(default_factory=set)
    description: Set[str] = field(default_factory=set)
    ncbi_gene: Set[str] = field(default_factory=set)
    ncbi_protein: Set[str] = field(default_factory=set)
    genbank_accession: Set[str] = field(default_factory=set)
    refseq_accession: Set[str] = field(default_factory=set)
    uniprot_accession: Set[str] = field(default_factory=set)
    sequence: Set[str] = field(default_factory=set)
    strand: Set[str] = field(default_factory=set)
    position_left: Set[int] = field(default_factory=set)
    position_right: Set[int] = field(default_factory=set)
    annotation_score: int = 0

    def to_df(self):
        pass


@dataclass
class EffectorDTO:
    name: Set[str] = field(default_factory=set)
    synonyms: Set[str] = field(default_factory=set)
    mechanism: Set[str] = field(default_factory=set)
    kegg_compounds: Set[str] = field(default_factory=set)

    def to_df(self):
        pass


@dataclass
class OrganismDTO:
    name: Set[str] = field(default_factory=set)
    species: Set[str] = field(default_factory=set)
    strain: Set[str] = field(default_factory=set)
    family: Set[str] = field(default_factory=set)
    phylum: Set[str] = field(default_factory=set)
    ncbi_taxonomy: Set[str] = field(default_factory=set)
    refseq_accession: Set[str] = field(default_factory=set)

    def to_df(self):
        pass


@dataclass
class PathwayDTO:
    name: Set[str] = field(default_factory=set)
    synonyms: Set[str] = field(default_factory=set)
    kegg_pathways: Set[str] = field(default_factory=set)

    def to_df(self):
        pass


@dataclass
class PublicationDTO:
    doi: Set[str] = field(default_factory=set)
    pmid: Set[str] = field(default_factory=set)
    title: Set[str] = field(default_factory=set)
    author: Set[str] = field(default_factory=set)
    year: Set[str] = field(default_factory=set)

    def to_df(self):
        pass
