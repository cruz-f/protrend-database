from typing import List, Type

from protrend.bioapis.publication import PubMedPublication
from protrend.transform.dto import PublicationDTO
from protrend.utils.miscellaneous import args_length


def _fetch_publications(identifiers: List[str],
                        cls: Type[PubMedPublication]) -> List[PubMedPublication]:
    compounds = []

    for identifier in identifiers:
        compound = cls(identifier=identifier)
        compound.fetch()
        compounds.append(compound)

    return compounds


def _annotate_publication(pubmed_publication: PubMedPublication, publication_dto: PublicationDTO):

    if pubmed_publication.identifier:
        publication_dto.pmid.add(pubmed_publication.pmid)
        publication_dto.doi.add(pubmed_publication.doi)
        publication_dto.title.add(pubmed_publication.title)
        publication_dto.author.add(pubmed_publication.author)
        publication_dto.year.add(pubmed_publication.year)


def annotate_publications(dtos: List[PublicationDTO],
                          identifiers: List[str] = None) -> PublicationDTO:
    """
    A common method to annotate a given effector with relevant information from NCBI PubMed.

    A given publication is annotated as follows:

        - 2º Step:
            - query pubmed identifiers (pmid) to NCBI PubMed database
            - retrieve publication relevant attributes

    :type dtos: List[PublicationDTO]
    :type identifiers: List[str]

    :rtype: List[PublicationDTO]

    :param dtos: list of PublicationDTO. Each publication DTO will be annotated with information from NCBI PubMed
    :param identifiers: list of identifiers to assist with the annotation

    :return: list of annotated PublicationDTO. This function returns the same list object for convenience
    """
    dtos_size = len(dtos)

    size = args_length(identifiers)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    pubmed_publications = _fetch_publications(identifiers=identifiers,
                                              cls=PubMedPublication)

    for publication_dto, pubmed_publication in zip(dtos, pubmed_publications):
        _annotate_publication(pubmed_publication, publication_dto)

    return dtos
