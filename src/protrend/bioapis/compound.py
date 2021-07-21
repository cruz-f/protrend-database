import os
from typing import List

import whoosh.index as w_index

from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.kegg import search_kegg_list


class Compound(BioAPI):

    def __init__(self, identifier: str = '', name: str = ''):

        super().__init__(identifier)

        self._name = name
        self._kegg_compounds = []

    @property
    def identifier(self) -> str:
        if self.kegg_compounds:
            return self.kegg_compounds[0]

        return self.identifier

    @property
    def kegg_compounds(self) -> List[str]:
        return self._kegg_compounds

    @property
    def name(self) -> str:
        return self._name

    def fetch(self, index: w_index.FileIndex = None):

        if index is None:

            cdw = os.getcwd()
            index_dir = os.path.join(cdw, 'compound_index')

            if os.path.exists(index_dir):
                index = w_index.open_dir(index_dir)

        if index is None:
            raise RuntimeError('Could not found a valid index for KEGG Compound')

        results_set = search_kegg_list(index=index, query=self.name, db='compound')

        if results_set:
            self._kegg_compounds.extend(list(results_set))