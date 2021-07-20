from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.entrez import entrez_summary, entrez_search


class Organism(BioAPI):

    def __init__(self, taxonomy: str = None, name: str = None):

        super().__init__()

        self._taxonomy = taxonomy
        self._name = name
        self._assembly_record = {}

    @property
    def taxonomy(self):
        return self.record.get('Id', self._taxonomy)

    @property
    def name(self):

        return self.record.get('ScientificName', self._name)

    @property
    def species(self):

        names = self.record.get('ScientificName', '').split()

        if names:
            return ' '.join(names[:2])

    @property
    def strain(self):

        names = self.record.get('ScientificName', '').split()

        if names:
            return ' '.join(names[2:])

    @property
    def assembly_record(self):
        return self._assembly_record

    @assembly_record.setter
    def assembly_record(self, value):
        if value:
            self._assembly_record = value

    @property
    def refseq(self):
        return self.assembly_record.get('AssemblyAccession')

    @property
    def refseq_ftp(self):
        return self.assembly_record.get('FtpPath_RefSeq')

    @property
    def genbank(self):

        return self.assembly_record.get('Synonym', {}).get('Genbank')

    @property
    def genbank_ftp(self):
        return self.assembly_record.get('FtpPath_GenBank')

    @property
    def assembly(self):
        if hasattr(self.assembly_record, 'attributes'):
            return self.assembly_record.attributes.get('uid')

    @property
    def assembly_name(self):
        return self.assembly_record.get('AssemblyName')

    def get_taxonomy_record(self):

        if self._taxonomy:
            self.record = entrez_summary(db='taxonomy', identifier=self._taxonomy)

        else:

            term = f'{self._name}[Scientific Name]'

            search_record = entrez_search(db='taxonomy', term=term)

            id_list = search_record.get('IdList', [])

            if id_list:
                identifier = id_list[0]
                self.record = entrez_summary(db='taxonomy', identifier=identifier)

    def fetch(self):

        self.get_taxonomy_record()

        assembly_term = f'txid{self.taxonomy}[Organism:noexp] AND "reference genome"[refseq category]'

        assembly_search_record = entrez_search(db='assembly', term=assembly_term)

        id_list = assembly_search_record.get('IdList', [])

        # broad search. However, it can retrieve multiple assemblies for the same organism. In this case, the first one
        # is selected.
        if not id_list:

            assembly_term = f'txid{self.taxonomy}[Organism:noexp]'

            assembly_search_record = entrez_search(db='assembly', term=assembly_term)

            id_list = assembly_search_record.get('IdList', [])

        if id_list:
            identifier = id_list[0]
            self.assembly_record = entrez_summary(db='assembly', identifier=identifier)
