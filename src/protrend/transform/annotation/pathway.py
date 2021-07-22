from typing import List, Type

import whoosh.index as w_index

from protrend.bioapis.pathway import KEGGPathway
from protrend.bioapis.kegg import indexing_kegg_list, fetch_kegg_list
from protrend.bioapis.utils import BIO_APIS_DIR
from protrend.transform.dto import PathwayDTO
from protrend.utils.miscellaneous import args_length


def _fetch_pathways(names: List[str],
                    cls: Type[KEGGPathway],
                    index: w_index.FileIndex = None) -> List[KEGGPathway]:
    compounds = []

    for name in names:
        compound = cls(name=name)
        compound.fetch()
        compounds.append(compound)

    return compounds


def _annotate_pathway(kegg_pathway: KEGGPathway, pathway_dto: PathwayDTO):
    pathway_dto.name.add(kegg_pathway.name)

    if kegg_pathway.identifier:
        pathway_dto.kegg_pathways.update(kegg_pathway.kegg_identifiers)
        pathway_dto.synonyms.update(kegg_pathway.kegg_names)


def annotate_pathways(dtos: List[PathwayDTO],
                      names: List[str] = None) -> PathwayDTO:
    """
    A common method to annotate a given effector with relevant information from KEGG Pathway.

    A given effector is annotated as follows:

        - 1º Step:
            - obtain KEGG list for pathway database
            - text indexing for all KEGG pathway names and identifiers

        - 2º Step:
            - query names to KEGG Pathway index
            - retrieve KEGG Pathway identifiers

    :type dtos: List[PathwayDTO]
    :type names: List[str]

    :rtype: List[PathwayDTO]

    :param dtos: list of PathwayDTO. Each pathway DTO will be annotated with information from KEGG
    :param names: list of names to assist with the annotation

    :return: list of annotated PathwayDTO. This function returns the same list object for convenience
    """
    dtos_size = len(dtos)

    size = args_length(names)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    db = 'pathway'
    cache = True
    try:
        index = indexing_kegg_list(db=db, cache=cache)

    except ValueError:

        print(f'Index for kegg database {db} was not found in path {BIO_APIS_DIR}. '
              f'Downloading kegg list and indexing ...')

        df_kegg_list = fetch_kegg_list(db=db, cache=cache)
        index = indexing_kegg_list(db=db, df_kegg_list=df_kegg_list)

    kegg_pathways = _fetch_pathways(names=names,
                                    cls=KEGGPathway,
                                    index=index)

    for pathway_dto, kegg_pathway in zip(dtos, kegg_pathways):
        _annotate_pathway(kegg_pathway, pathway_dto)

    return dtos
