import os
from typing import List

import whoosh.index as w_index

from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.kegg import search_kegg_list, KEGG_PATH
from protrend.utils.miscellaneous import is_null


class KEGGPathway(BioAPI):

    def __init__(self, identifier: str = '', name: str = ''):

        super().__init__(identifier)

        if is_null(name):
            name = ''

        self._name = name
        self._kegg_identifiers = []
        self._kegg_names = []

    @property
    def identifier(self) -> str:
        if self.kegg_identifiers:
            return self.kegg_identifiers[0]

        return self._identifier

    @property
    def name(self) -> str:
        return self._name

    @property
    def kegg_identifiers(self) -> List[str]:
        return self._kegg_identifiers

    @property
    def kegg_names(self) -> List[str]:
        return self._kegg_names

    def fetch(self, index: w_index.FileIndex = None):

        if index is None:

            index_dir = os.path.join(KEGG_PATH, 'pathway_index')

            if os.path.exists(index_dir):
                index = w_index.open_dir(index_dir)

        if index is None:
            raise RuntimeError('Could not found a valid index for KEGG Pathway')

        if self.name:
            identifiers, names = search_kegg_list(index=index, query=self.name, db='pathway')

            if identifiers:
                self._kegg_identifiers.extend(list(identifiers))

            if names:
                self._kegg_names.extend(list(names))