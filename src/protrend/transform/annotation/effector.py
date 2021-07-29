from typing import List, Type

import whoosh.index as w_index

from protrend.bioapis.compound import KEGGCompound
from protrend.bioapis.kegg import indexing_kegg_list, fetch_kegg_list
from protrend.bioapis.utils import BIO_APIS_DIR
from protrend.transform.dto import EffectorDTO
from protrend.utils.miscellaneous import args_length


def _fetch_compounds(names: List[str],
                     cls: Type[KEGGCompound],
                     index: w_index.FileIndex = None) -> List[KEGGCompound]:
    compounds = []

    for name in names:
        compound = cls(name=name)
        compound.fetch()
        compounds.append(compound)

    return compounds


def _annotate_compound(kegg_compound: KEGGCompound, effector_dto: EffectorDTO):

    effector_dto.name.append(kegg_compound.name)

    if kegg_compound.identifier:
        effector_dto.kegg_compounds.extend(kegg_compound.kegg_identifiers)
        effector_dto.synonyms.extend(kegg_compound.kegg_names)


def annotate_effectors(dtos: List[EffectorDTO],
                       names: List[str] = None) -> EffectorDTO:
    """
    A common method to annotate a given effector with relevant information from KEGG Compound.

    A given effector is annotated as follows:

        - 1º Step:
            - obtain KEGG list for compound database
            - text indexing for all KEGG compound names and identifiers

        - 2º Step:
            - query names to KEGG compound index
            - retrieve KEGG Compound identifiers

    :type dtos: List[EffectorDTO]
    :type names: List[str]

    :rtype: List[EffectorDTO]

    :param dtos: list of EffectorDTO. Each effector DTO will be annotated with information from KEGG
    :param names: list of names to assist with the annotation

    :return: list of annotated EffectorDTO. This function returns the same list object for convenience
    """
    dtos_size = len(dtos)

    size = args_length(names)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    db = 'compound'
    cache = True
    try:
        index = indexing_kegg_list(db=db, cache=cache)

    except ValueError:

        print(f'Index for kegg database {db} was not found in path {BIO_APIS_DIR}. '
              f'Downloading kegg list and indexing ...')

        df_kegg_list = fetch_kegg_list(db=db, cache=cache)
        index = indexing_kegg_list(db=db, df_kegg_list=df_kegg_list)

    kegg_compounds = _fetch_compounds(names=names,
                                      cls=KEGGCompound,
                                      index=index)

    for effector_dto, kegg_compound in zip(dtos, kegg_compounds):

        _annotate_compound(kegg_compound, effector_dto)

    return dtos
