from typing import List, Type, TYPE_CHECKING

from tqdm import tqdm

from protrend.bioapis import NCBITaxonomyOrganism
from protrend.log import ProtrendLogger
from protrend.utils.miscellaneous import args_length, scale_arg


if TYPE_CHECKING:
    from .dto import OrganismDTO


def _fetch_organisms(identifiers: List[str],
                     names: List[str],
                     cls: Type[NCBITaxonomyOrganism]) -> List[NCBITaxonomyOrganism]:
    organisms = []

    for identifier, name in tqdm(zip(identifiers, names), desc='organism', total=len(identifiers)):
        organism = cls(identifier=identifier, name=name)
        organism.fetch()
        organisms.append(organism)

    return organisms


def _annotate_organism(ncbi_taxonomy: NCBITaxonomyOrganism, organism_dto: 'OrganismDTO'):
    if ncbi_taxonomy.identifier:
        organism_dto.name.append(ncbi_taxonomy.name)
        organism_dto.species.append(ncbi_taxonomy.species)
        organism_dto.strain.append(ncbi_taxonomy.strain)
        organism_dto.ncbi_taxonomy.append(ncbi_taxonomy.taxonomy)
        organism_dto.refseq_accession.append(ncbi_taxonomy.refseq)
        organism_dto.refseq_ftp.append(ncbi_taxonomy.refseq_ftp)
        organism_dto.genbank_accession.append(ncbi_taxonomy.genbank)
        organism_dto.genbank_ftp.append(ncbi_taxonomy.genbank_ftp)
        organism_dto.ncbi_assembly.append(ncbi_taxonomy.assembly)
        organism_dto.assembly_accession.append(ncbi_taxonomy.assembly_name)


def annotate_organisms(dtos: List['OrganismDTO'],
                       identifiers: List[str] = None,
                       names: List[str] = None) -> List['OrganismDTO']:
    """
    A common method to annotate a given organism with relevant information from NCBI Taxonomy.

    A given organism is annotated as follows:

        - 1º Step:
            - query organism names or NCBI taxonomy identifiers to NCBI Taxonomy
            - retrieve the NCBI taxonomy record for the organism
            - query NCBI Assembly database using the refseq genome assembly whenever possible
            - retrieve the NCBI assembly record for the organism
            - retrieve organism relevant attributes from taxonomy and assembly records

    :type dtos: List['OrganismDTO']
    :type identifiers: List[str]
    :type names: List[str]

    :param dtos: list of 'OrganismDTO'. Each publication DTO will be annotated with information from NCBI PubMed
    :param identifiers: list of identifiers to assist with the annotation
    :param names: list of names to assist with the annotation

    :return: list of annotated 'OrganismDTO'. This function returns the same list object for convenience
    """
    if not dtos:
        return dtos

    dtos_size = len(dtos)

    size = args_length(identifiers, names)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    identifiers = scale_arg(identifiers, size)
    names = scale_arg(names, size)

    ProtrendLogger.log.info(f'Starting fetch {len(names)} organisms to cls {NCBITaxonomyOrganism.__name__}')
    ncbi_taxonomys = _fetch_organisms(identifiers=identifiers,
                                      names=names,
                                      cls=NCBITaxonomyOrganism)
    ProtrendLogger.log.info(f'Finishing fetch organisms')

    for organism_dto, ncbi_taxonomy in zip(dtos, ncbi_taxonomys):
        _annotate_organism(ncbi_taxonomy, organism_dto)

    return dtos
