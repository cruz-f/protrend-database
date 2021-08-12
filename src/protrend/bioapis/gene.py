from typing import List

from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.entrez import entrez_summary, entrez_search


class NCBIGene(BioAPI):

    def __init__(self,
                 identifier: str = '',
                 taxonomy: int = 0,
                 locus_tag: str = '',
                 name: str = ''):

        super().__init__(identifier)

        self._taxonomy = taxonomy
        self._locus_tag = locus_tag
        self._name = name

    @property
    def identifier(self) -> str:
        if hasattr(self.record, 'attributes'):
            return self.record.attributes.get('uid', self._identifier)

        return self._identifier

    @property
    def taxonomy(self) -> int:
        return self.record.get('Organism', {}).get('TaxID', self._taxonomy)

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
        return self.record.get('OtherAliases', '').split(', ')

    @property
    def refseq_accession(self) -> str:

        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:
            genomic_info = genomic_infos[0]

            return genomic_info.get('ChrAccVer', '')

        return ''

    @property
    def position_left(self) -> int:
        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:
            genomic_info = genomic_infos[0]

            left = genomic_info.get('ChrStart', None)

            if left is not None:
                return int(left)

    @property
    def position_right(self) -> int:
        genomic_infos = self.record.get('GenomicInfo')

        if genomic_infos:
            genomic_info = genomic_infos[0]

            right = genomic_info.get('ChrStop', None)

            if right is not None:
                return int(right)

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
                search_record, _ = entrez_search(db='gene', term=term)

                id_list = search_record.get('IdList', [])

                if id_list:
                    identifier = id_list[0]

        if identifier:
            record, _ = entrez_summary(db='gene', identifier=identifier)
            self.record = record
