from typing import List, Type, TYPE_CHECKING

import whoosh.index as w_index
from tqdm import tqdm

from protrend.bioapis import KEGGPathway, indexing_kegg_list, fetch_kegg_list
from protrend.log import ProtrendLogger
from protrend.utils.miscellaneous import args_length
from ..utils import Settings

if TYPE_CHECKING:
    from .dto import PathwayDTO


def _fetch_pathways(names: List[str],
                    cls: Type[KEGGPathway],
                    index: w_index.FileIndex = None) -> List[KEGGPathway]:
    pathways = []

    for name in tqdm(names, desc='pathway'):
        pathway = cls(name=name)
        pathway.fetch(index)
        pathways.append(pathway)

    return pathways


def _annotate_pathway(kegg_pathway: KEGGPathway, pathway_dto: 'PathwayDTO'):
    pathway_dto.name.append(kegg_pathway.name)

    if kegg_pathway.identifier:
        pathway_dto.kegg_pathways.extend(kegg_pathway.kegg_identifiers)
        pathway_dto.synonyms.extend(kegg_pathway.kegg_names)


def annotate_pathways(dtos: List['PathwayDTO'],
                      names: List[str] = None) -> List['PathwayDTO']:
    """
    A common method to annotate a given pathway with relevant information from KEGG Pathway.

    A given pathway is annotated as follows:

        - 1º Step:
            - obtain KEGG list for pathway database
            - text indexing for all KEGG pathway names and identifiers

        - 2º Step:
            - query names to KEGG Pathway index
            - retrieve KEGG Pathway identifiers

    :type dtos: List['PathwayDTO']
    :type names: List[str]

    :rtype: List['PathwayDTO']

    :param dtos: list of 'PathwayDTO'. Each pathway DTO will be annotated with information from KEGG
    :param names: list of names to assist with the annotation

    :return: list of annotated 'PathwayDTO'. This function returns the same list object for convenience
    """
    if not dtos:
        return dtos

    dtos_size = len(dtos)

    size = args_length(names)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    db = 'pathway'
    # noinspection DuplicatedCode
    cache = True
    try:
        index = indexing_kegg_list(db=db, cache=cache)

    except ValueError:

        ProtrendLogger.log.debug(f'Index for kegg database {db} was not found in path {Settings.kegg}. '
                                 f'Downloading kegg list and indexing ...')

        df_kegg_list = fetch_kegg_list(db=db, cache=cache)
        index = indexing_kegg_list(db=db, df_kegg_list=df_kegg_list)

    ProtrendLogger.log.info(f'Starting fetch {len(names)} pathways to cls {KEGGPathway.__name__}')
    kegg_pathways = _fetch_pathways(names=names, cls=KEGGPathway, index=index)
    ProtrendLogger.log.info(f'Finishing fetch pathways')

    for pathway_dto, kegg_pathway in zip(dtos, kegg_pathways):
        _annotate_pathway(kegg_pathway, pathway_dto)

    return dtos
