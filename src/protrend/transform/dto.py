from dataclasses import dataclass, field

from protrend.utils.unique_list import UniqueList


@dataclass
class GeneDTO:
    locus_tag: UniqueList[str] = field(default_factory=UniqueList)
    name: UniqueList[str] = field(default_factory=UniqueList)
    synonyms: UniqueList[str] = field(default_factory=UniqueList)
    function: UniqueList[str] = field(default_factory=UniqueList)
    description: UniqueList[str] = field(default_factory=UniqueList)
    ncbi_gene: UniqueList[str] = field(default_factory=UniqueList)
    ncbi_protein: UniqueList[str] = field(default_factory=UniqueList)
    genbank_accession: UniqueList[str] = field(default_factory=UniqueList)
    refseq_accession: UniqueList[str] = field(default_factory=UniqueList)
    uniprot_accession: UniqueList[str] = field(default_factory=UniqueList)
    sequence: UniqueList[str] = field(default_factory=UniqueList)
    strand: UniqueList[str] = field(default_factory=UniqueList)
    position_left: UniqueList[int] = field(default_factory=UniqueList)
    position_right: UniqueList[int] = field(default_factory=UniqueList)
    annotation_score: int = 0

    def to_df(self):
        pass


@dataclass
class EffectorDTO:
    name: UniqueList[str] = field(default_factory=UniqueList)
    synonyms: UniqueList[str] = field(default_factory=UniqueList)
    mechanism: UniqueList[str] = field(default_factory=UniqueList)
    kegg_compounds: UniqueList[str] = field(default_factory=UniqueList)

    def to_df(self):
        pass


@dataclass
class OrganismDTO:
    name: UniqueList[str] = field(default_factory=UniqueList)
    species: UniqueList[str] = field(default_factory=UniqueList)
    strain: UniqueList[str] = field(default_factory=UniqueList)
    family: UniqueList[str] = field(default_factory=UniqueList)
    phylum: UniqueList[str] = field(default_factory=UniqueList)
    ncbi_taxonomy: UniqueList[str] = field(default_factory=UniqueList)
    refseq_accession: UniqueList[str] = field(default_factory=UniqueList)
    refseq_ftp: UniqueList[str] = field(default_factory=UniqueList)
    genbank_accession: UniqueList[str] = field(default_factory=UniqueList)
    genbank_ftp: UniqueList[str] = field(default_factory=UniqueList)
    ncbi_assembly: UniqueList[str] = field(default_factory=UniqueList)
    assembly_accession: UniqueList[str] = field(default_factory=UniqueList)

    def to_df(self):
        pass


@dataclass
class PathwayDTO:
    name: UniqueList[str] = field(default_factory=UniqueList)
    synonyms: UniqueList[str] = field(default_factory=UniqueList)
    kegg_pathways: UniqueList[str] = field(default_factory=UniqueList)

    def to_df(self):
        pass


@dataclass
class PublicationDTO:
    pmid: UniqueList[str] = field(default_factory=UniqueList)
    doi: UniqueList[str] = field(default_factory=UniqueList)
    title: UniqueList[str] = field(default_factory=UniqueList)
    author: UniqueList[str] = field(default_factory=UniqueList)
    year: UniqueList[str] = field(default_factory=UniqueList)

    def to_df(self):
        pass
