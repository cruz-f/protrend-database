from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.entrez import entrez_summary, entrez_search
from protrend.utils.miscellaneous import is_null


class NCBITaxonomyOrganism(BioAPI):

    def __init__(self, identifier: str = '', name: str = ''):

        super().__init__(identifier)

        if is_null(name):
            name = ''

        name = str(name)

        self._name = name
        self._assembly_record = {}

    @property
    def identifier(self) -> str:
        return str(self.record.get('Id', self._identifier))

    @property
    def taxonomy(self) -> str:
        return str(self.identifier)

    @property
    def name(self) -> str:

        return self.record.get('ScientificName', self._name)

    @property
    def species(self) -> str:

        names = self.record.get('ScientificName', '').split()

        if names:
            return ' '.join(names[:2])

        return ''

    @property
    def strain(self) -> str:

        names = self.record.get('ScientificName', '').split()

        if names:
            return ' '.join(names[2:])

        return ''

    @property
    def assembly_record(self) -> dict:
        return self._assembly_record

    @assembly_record.setter
    def assembly_record(self, value):
        if value:
            self._assembly_record = value

    @property
    def refseq(self) -> str:
        return self.assembly_record.get('AssemblyAccession', '')

    @property
    def refseq_ftp(self) -> str:
        return self.assembly_record.get('FtpPath_RefSeq', '')

    @property
    def genbank(self) -> str:

        return self.assembly_record.get('Synonym', {}).get('Genbank', '')

    @property
    def genbank_ftp(self) -> str:
        return self.assembly_record.get('FtpPath_GenBank', '')

    @property
    def assembly(self) -> str:
        if hasattr(self.assembly_record, 'attributes'):
            return self.assembly_record.attributes.get('uid', '')

        return ''

    @property
    def assembly_name(self) -> str:
        return self.assembly_record.get('AssemblyName', '')

    def get_taxonomy_record(self):

        if self._identifier:
            record = entrez_summary(db='taxonomy', identifier=self._identifier)
            self.record = record

        else:

            if self._name:
                term = f'{self._name}[Scientific Name]'

                search_record = entrez_search(db='taxonomy', term=term)

                id_list = search_record.get('IdList', [])

                if id_list:
                    identifier = id_list[0]
                    record = entrez_summary(db='taxonomy', identifier=identifier)
                    self.record = record

    def fetch(self):

        self.get_taxonomy_record()

        if self.taxonomy:
            assembly_term = f'txid{self.taxonomy}[Organism:noexp] AND "reference genome"[refseq category]'

            assembly_search_record = entrez_search(db='assembly', term=assembly_term)

            id_list = assembly_search_record.get('IdList', [])

            # broad search. However, it can retrieve multiple assemblies for the same organism.
            # In this case, the first one is selected.
            if not id_list:

                assembly_term = f'txid{self.taxonomy}[Organism:noexp]'

                assembly_search_record = entrez_search(db='assembly', term=assembly_term)

                id_list = assembly_search_record.get('IdList', [])

            if id_list:
                identifier = id_list[0]
                record = entrez_summary(db='assembly', identifier=identifier)
                self.assembly_record = record
