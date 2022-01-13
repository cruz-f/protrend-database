from typing import List, Type, TYPE_CHECKING

import whoosh.index as w_index
from tqdm import tqdm

from protrend.bioapis import KEGGCompound
from protrend.bioapis import indexing_kegg_list, fetch_kegg_list
from protrend.log import ProtrendLogger
from protrend.utils.miscellaneous import args_length
from ..utils import Settings

if TYPE_CHECKING:
    from .dto import EffectorDTO


def _fetch_compounds(names: List[str],
                     cls: Type[KEGGCompound],
                     index: w_index.FileIndex = None) -> List[KEGGCompound]:
    compounds = []

    for name in tqdm(names, desc='compound'):
        compound = cls(name=name)
        compound.fetch(index)
        compounds.append(compound)

    return compounds


def _annotate_compound(kegg_compound: KEGGCompound, effector_dto: 'EffectorDTO'):
    effector_dto.name.append(kegg_compound.name)

    if kegg_compound.identifier:
        effector_dto.kegg_compounds.extend(kegg_compound.kegg_identifiers)
        effector_dto.synonyms.extend(kegg_compound.kegg_names)


def annotate_effectors(dtos: List['EffectorDTO'],
                       names: List[str] = None) -> List['EffectorDTO']:
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
    if not dtos:
        return dtos

    dtos_size = len(dtos)

    size = args_length(names)

    if size != dtos_size:
        raise ValueError(f'Invalid inputs for dto list size of {dtos_size} and args size of {size}')

    db = 'compound'
    # noinspection DuplicatedCode
    cache = True
    try:
        index = indexing_kegg_list(db=db, cache=cache)

    except ValueError:

        ProtrendLogger.log.debug(f'Index for kegg database {db} was not found in path {Settings.kegg}. '
                                 f'Downloading kegg list and indexing ...')

        df_kegg_list = fetch_kegg_list(db=db, cache=cache)
        index = indexing_kegg_list(db=db, df_kegg_list=df_kegg_list)

    ProtrendLogger.log.info(f'Starting fetch {len(names)} compounds to cls {KEGGCompound.__name__}')
    kegg_compounds = _fetch_compounds(names=names, cls=KEGGCompound, index=index)
    ProtrendLogger.log.info(f'Finishing fetch compounds')

    for effector_dto, kegg_compound in zip(dtos, kegg_compounds):
        _annotate_compound(kegg_compound, effector_dto)

    return dtos
