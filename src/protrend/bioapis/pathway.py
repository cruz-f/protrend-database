import os

import whoosh.index as w_index

from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.kegg import search_kegg_list


class Pathway(BioAPI):

    def __init__(self, identifier: str = '', name: str = ''):

        super().__init__(identifier)

        self._name = name
        self._kegg_pathways = []

    @property
    def identifier(self):
        if self.kegg_pathways:
            return self.kegg_pathways[0]

        return self.identifier

    @property
    def kegg_pathways(self):
        return self._kegg_pathways

    @property
    def name(self):
        return self._name

    def fetch(self, index: w_index.FileIndex = None):

        if index is None:

            cdw = os.getcwd()
            index_dir = os.path.join(cdw, 'pathway_index')

            if os.path.exists(index_dir):
                index = w_index.open_dir(index_dir)

        if index is None:
            raise RuntimeError('Could not found a valid index for KEGG Pathway')

        results_set = search_kegg_list(index=index, query=self.name, db='pathway')

        if results_set:
            self._kegg_pathways.extend(list(results_set))