from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.entrez import entrez_summary, entrez_search


class NCBIGene(BioAPI):

    def __init__(self,
                 gene: str = None,
                 taxonomy: str = None,
                 locus_tag: str = None,
                 name: str = None):

        super().__init__()

        self._gene = gene
        self._taxonomy = taxonomy
        self._locus_tag = locus_tag
        self._name = name

    @property
    def gene(self):
        if hasattr(self.record, 'attributes'):
            return self.record.attributes.get('uid', self._gene)

        return self._gene

    @property
    def taxonomy(self):
        return self.record.get('Organism', {}).get('TaxID', self._taxonomy)

    @property
    def locus_tag(self):

        if self.synonyms:
            return self.synonyms[0]

        return self._locus_tag

    @property
    def name(self):
        return self.record.get('Name', self._name)

    @property
    def function(self):
        return self.record.get('Description')

    @property
    def synonyms(self):
        return self.record.get('OtherAliases', '').split(', ')

    @property
    def refseq_accession(self):

        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:

            genomic_info = genomic_infos[0]

            return genomic_info.get('ChrAccVer')

    @property
    def position_left(self):
        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:
            genomic_info = genomic_infos[0]

            left = genomic_info.get('ChrStart', None)

            if left is not None:
                return int(left)

    @property
    def position_right(self):
        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:
            genomic_info = genomic_infos[0]

            right = genomic_info.get('ChrStop', None)

            if right is not None:
                return int(right)

    def build_term(self):

        if self._locus_tag and self._taxonomy:

            term = f'{self._locus_tag} AND txid{self._taxonomy}[Organism:noexp]'

        elif self._locus_tag and not self._taxonomy:

            term = f'{self._locus_tag}'

        elif self._name and self._taxonomy:

            term = f'{self._name} AND txid{self._taxonomy}[Organism:noexp]'

        else:
            term = ''

        return term

    def fetch(self):

        if self._gene:
            identifier = self._gene

        else:
            identifier = None

            term = self.build_term()

            if term:
                search_record = entrez_search(db='gene', term=term)

                id_list = search_record.get('IdList', [])

                if id_list:
                    identifier = id_list[0]

        if identifier:
            self.record = entrez_summary(db='gene', identifier=identifier)
