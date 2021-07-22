from typing import List, Type

from protrend.bioapis.organism import NCBITaxonomyOrganism
from protrend.transform.dto import OrganismDTO
from protrend.utils.miscellaneous import args_length


def _fetch_organisms(identifiers: List[str],
                     names: List[str],
                     cls: Type[NCBITaxonomyOrganism]) -> List[NCBITaxonomyOrganism]:
    organisms = []

    for identifier, name in zip(identifiers, names):
        organism = cls(identifier=identifier, name=name)
        organism.fetch()
        organisms.append(organism)

    return organisms


def _annotate_organism(ncbi_taxonomy: NCBITaxonomyOrganism, organism_dto: OrganismDTO):
    if ncbi_taxonomy.identifier:
        organism_dto.name.add(ncbi_taxonomy.name)
        organism_dto.species.add(ncbi_taxonomy.species)
        organism_dto.strain.add(ncbi_taxonomy.strain)
        organism_dto.ncbi_taxonomy.add(ncbi_taxonomy.taxonomy)
        organism_dto.refseq_accession.add(ncbi_taxonomy.refseq)
        organism_dto.refseq_ftp.add(ncbi_taxonomy.refseq_ftp)
        organism_dto.genbank_accession.add(ncbi_taxonomy.genbank)
        organism_dto.genbank_ftp.add(ncbi_taxonomy.genbank_ftp)
        organism_dto.ncbi_assembly.add(ncbi_taxonomy.assembly)
        organism_dto.assembly_accession.add(ncbi_taxonomy.assembly_name)


def annotate_organisms(dtos: List[OrganismDTO],
                       identifiers: List[str] = None,
                       names: List[str] = None) -> OrganismDTO:
    """
    A common method to annotate a given organism with relevant information from NCBI Taxonomy.

    A given organism is annotated as follows:

        - 1º Step:
            - query organism names or NCBI taxonomy identifiers to NCBI Taxonomy
            - retrieve the NCBI taxonomy record for the organism
            - query NCBI Assembly database using the refseq genome assembly whenever possible
            - retrieve the NCBI assembly record for the organism
            - retrieve organism relevant attributes from taxonomy and assembly records

    :type dtos: List[OrganismDTO]
    :type identifiers: List[str]
    :type names: List[str]

    :rtype: List[OrganismDTO]

    :param dtos: list of OrganismDTO. Each publication DTO will be annotated with information from NCBI PubMed
    :param identifiers: list of identifiers to assist with the annotation
    :param names: list of names to assist with the annotation

    :return: list of annotated OrganismDTO. This function returns the same list object for convenience
    """
    dtos_size = len(dtos)

    size = args_length(identifiers, names)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    ncbi_taxonomys = _fetch_organisms(identifiers=identifiers,
                                      names=names,
                                      cls=NCBITaxonomyOrganism)

    for organism_dto, ncbi_taxonomy in zip(dtos, ncbi_taxonomys):
        _annotate_organism(ncbi_taxonomy, organism_dto)

    return dtos
