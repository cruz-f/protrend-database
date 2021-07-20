from protrend.bioapis.bioapi import BioAPI


class Metabolite(BioAPI):

    def __init__(self, name: str = None):

        super().__init__()

        self._name = name

    @property
    def kegg_metabolite(self):
        return 1