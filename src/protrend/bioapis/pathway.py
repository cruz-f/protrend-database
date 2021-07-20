from protrend.bioapis.bioapi import BioAPI


class Pathway(BioAPI):

    def __init__(self, name: str = None):

        super().__init__()

        self._name = name

    @property
    def kegg_pathway(self):
        return 1