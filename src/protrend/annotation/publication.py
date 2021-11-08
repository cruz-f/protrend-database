from typing import List, Type, TYPE_CHECKING

from protrend.bioapis import PubMedPublication
from protrend.log import ProtrendLogger
from protrend.utils.miscellaneous import args_length


if TYPE_CHECKING:
    from .dto import PublicationDTO


def _fetch_publications(identifiers: List[str],
                        cls: Type[PubMedPublication]) -> List[PubMedPublication]:
    publications = []

    for identifier in identifiers:
        publication = cls(identifier=identifier)
        publication.fetch()
        publications.append(publication)

    return publications


def _annotate_publication(pubmed_publication: PubMedPublication, publication_dto: 'PublicationDTO'):
    if pubmed_publication.identifier:
        publication_dto.pmid.append(pubmed_publication.pmid)
        publication_dto.doi.append(pubmed_publication.doi)
        publication_dto.title.append(pubmed_publication.title)
        publication_dto.author.append(pubmed_publication.author)
        publication_dto.year.append(pubmed_publication.year)


def annotate_publications(dtos: List['PublicationDTO'],
                          identifiers: List[str] = None) -> List['PublicationDTO']:
    """
    A common method to annotate a given publication with relevant information from NCBI PubMed.

    A given publication is annotated as follows:

        - 1º Step:
            - query pubmed identifiers (pmid) to NCBI PubMed database
            - retrieve publication relevant attributes

    :type dtos: List['PublicationDTO']
    :type identifiers: List[str]

    :rtype: List['PublicationDTO']

    :param dtos: list of 'PublicationDTO'. Each publication DTO will be annotated with information from NCBI PubMed
    :param identifiers: list of identifiers to assist with the annotation

    :return: list of annotated 'PublicationDTO'. This function returns the same list object for convenience
    """
    dtos_size = len(dtos)

    size = args_length(identifiers)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    ProtrendLogger.log.info(f'Starting fetch {len(identifiers)} publications to cls {PubMedPublication.__name__}')
    pubmed_publications = _fetch_publications(identifiers=identifiers,
                                              cls=PubMedPublication)
    ProtrendLogger.log.info(f'Finishing fetch publications')

    for publication_dto, pubmed_publication in zip(dtos, pubmed_publications):
        _annotate_publication(pubmed_publication, publication_dto)

    return dtos
