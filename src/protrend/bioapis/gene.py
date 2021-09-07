from typing import List

from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.entrez import entrez_summary, entrez_search
from protrend.utils.miscellaneous import is_null


class NCBIGene(BioAPI):

    def __init__(self,
                 identifier: str = '',
                 taxonomy: str = '',
                 locus_tag: str = '',
                 name: str = ''):

        super().__init__(identifier)

        if is_null(taxonomy):
            taxonomy = ''

        taxonomy = str(taxonomy)

        if is_null(locus_tag):
            locus_tag = ''

        locus_tag = str(locus_tag)

        if is_null(name):
            name = ''

        name = str(name)

        self._taxonomy = taxonomy
        self._locus_tag = locus_tag
        self._name = name

    @property
    def identifier(self) -> str:
        return self.record.get('attributes', {}).get('uid', self._identifier)

    @property
    def taxonomy(self) -> str:
        return str(self.record.get('Organism', {}).get('TaxID', self._taxonomy))

    @property
    def locus_tag(self) -> str:

        if self.synonyms:
            return self.synonyms[0]

        return self._locus_tag

    @property
    def name(self) -> str:
        return self.record.get('Name', self._name)

    @property
    def function(self) -> str:
        return self.record.get('Description', '')

    @property
    def synonyms(self) -> List[str]:
        synonyms = self.record.get('OtherAliases', '').split(', ')

        if self._locus_tag:
            synonyms.append(self._locus_tag)

        if self._name:
            synonyms.append(self._name)

        if self.name:
            synonyms.append(self.name)

        return synonyms

    @property
    def refseq_accession(self) -> str:

        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:
            genomic_info = genomic_infos[0]

            return genomic_info.get('ChrAccVer', '')

        return ''

    @property
    def start(self) -> int:
        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:
            genomic_info = genomic_infos[0]

            left = genomic_info.get('ChrStart', None)

            if left is not None:
                return int(left)

    @property
    def stop(self) -> int:
        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:
            genomic_info = genomic_infos[0]

            right = genomic_info.get('ChrStop', None)

            if right is not None:
                return int(right)

    @property
    def strand(self) -> str:
        if self.start is None:
            return 'unknown'

        if self.stop is None:
            return 'unknown'

        if self.start < self.stop:
            return 'forward'
        return 'reverse'

    def build_term(self) -> str:

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

        if self._identifier:
            identifier = self._identifier

        else:
            identifier = None

            term = self.build_term()

            if term:
                search_record = entrez_search(db='gene', term=term)

                id_list = search_record.get('IdList', [])

                if id_list:
                    identifier = id_list[0]

        if identifier:
            record = entrez_summary(db='gene', identifier=identifier)
            self.record = record
