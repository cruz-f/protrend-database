from typing import List

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.entrez import entrez_summary, entrez_search, entrez_fetch
from protrend.bioapis.uniprot import fetch_uniprot_record, query_uniprot


class NCBIProtein(BioAPI):

    def __init__(self,
                 identifier: str = '',
                 refseq_accession: str = '',
                 genbank_accession: str = '',
                 taxonomy: int = 0,
                 locus_tag: str = '',
                 name: str = ''):

        super().__init__(identifier)

        self._refseq_accession = refseq_accession
        self._genbank_accession = genbank_accession
        self._taxonomy = taxonomy
        self._locus_tag = locus_tag
        self._name = name
        self._seq_record = SeqRecord(seq=Seq(''))

    @property
    def identifier(self) -> str:
        return self.record.get('Id', self._identifier)

    def is_refseq(self) -> bool:
        additional = self.record.get('Extra', '')
        extra_info = additional.split('|')

        for info in extra_info:
            if info.lower() == 'ref':
                return True

        return False

    @property
    def refseq_accession(self) -> str:

        if self.is_refseq():
            return self.record.get('AccessionVersion', self._refseq_accession)

        return ''

    def is_genbank(self) -> bool:
        additional = self.record.get('Extra', '')
        extra_info = additional.split('|')

        for info in extra_info:
            if info.lower() == 'gb':
                return True

        return False

    @property
    def genbank_accession(self) -> str:

        if self.is_genbank():
            return self.record.get('AccessionVersion', self._genbank_accession)

        return ''

    @property
    def taxonomy(self) -> int:
        return self.record.get('TaxId', self._taxonomy)

    @property
    def locus_tag(self) -> str:

        if self.seq_record:
            for feature in self.seq_record.features:

                if feature.type == 'CDS':
                    locus_tag = feature.qualifiers.get('locus_tag', self._locus_tag)

                    if locus_tag:
                        return locus_tag[0]

        return self._locus_tag

    @property
    def name(self) -> str:
        return self.record.get('Caption', self._name)

    @property
    def seq_record(self) -> SeqRecord:
        return self._seq_record

    @seq_record.setter
    def seq_record(self, value):
        if value:
            self._seq_record = value

    @property
    def synonyms(self) -> List[str]:

        synonyms = []

        if self.locus_tag:
            synonyms.append(self.locus_tag)

        if self.seq_record:
            for feature in self.seq_record.features:

                if feature.type == 'CDS':
                    old_locus_tag = feature.qualifiers.get('old_locus_tag')

                    if old_locus_tag:
                        synonyms.append(old_locus_tag[0])

                        break

        if self.name:
            synonyms.append(self.name)

        return synonyms

    @property
    def sequence(self) -> str:
        return str(self.seq_record.seq)

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

        elif self._genbank_accession:
            identifier = self._genbank_accession

        elif self._refseq_accession:
            identifier = self._refseq_accession

        else:
            identifier = None

            term = self.build_term()

            if term:
                search_record = entrez_search(db='protein', term=term)

                id_list = search_record.get('IdList', [])

                if id_list:
                    identifier = id_list[0]

        if identifier:
            self.record = entrez_summary(db='protein', identifier=identifier)
            self.seq_record = entrez_fetch(db='protein', identifier=identifier)


class UniProtProtein(BioAPI):

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
        return getattr(self.record, 'id', self._identifier)

    @property
    def taxonomy(self) -> int:
        dbxrefs = getattr(self.record, 'dbxrefs', None)

        if dbxrefs:
            for ref in dbxrefs:
                if 'NCBI Taxonomy:' in ref:
                    return ref.split(':')[1]

        return self._taxonomy

    @property
    def locus_tag(self) -> str:

        annotations = getattr(self.record, 'annotations', {})

        if annotations:
            loci = annotations.get('gene_name_ordered locus', [])

            for locus in loci:
                return locus

        return self._locus_tag

    @property
    def name(self) -> str:
        annotations = getattr(self.record, 'annotations', {})

        if annotations:
            return annotations.get('gene_name_primary', self._name)

        if self.locus_tag:
            return self.locus_tag

        return self._name

    @property
    def function(self) -> str:
        annotations = getattr(self.record, 'annotations', {})

        if annotations:
            funcs = annotations.get('recommendedName_fullName', [])

            for func in funcs:
                return func

        return ''

    @property
    def description(self) -> str:
        annotations = getattr(self.record, 'annotations', {})

        if annotations:
            funcs = annotations.get('comment_function', [])

            for func in funcs:
                return func

        return ''

    @property
    def synonyms(self) -> List[str]:

        synonyms = []

        if self.locus_tag:
            synonyms.append(self.locus_tag)

        if self.name:
            synonyms.append(self.name)

        annotations = getattr(self.record, 'annotations', {})

        if annotations:
            loci = annotations.get('gene_name_ordered locus', [])
            synonyms.extend(loci)

        return synonyms

    @property
    def sequence(self) -> str:
        return str(getattr(self.record, 'seq', ''))

    def parse_uniprot_query(self, query: pd.DataFrame):

        if self._taxonomy:
            tax_mask = query.loc[:, 'Organism ID'] == self._taxonomy

            query = query[tax_mask]

        if self._locus_tag:

            loci = query.loc[:, 'Gene names']
            loci_mask = []

            for value in loci:

                mask_value = False

                values = value.split()

                for text in values:
                    text = ''.join(letter for letter in text if letter.isalnum())

                    if self._locus_tag.lower() in text.lower() or text.lower() in self._locus_tag.lower():
                        mask_value = True
                        break

                loci_mask.append(mask_value)

            accessions = query.loc[loci_mask, 'Entry']

            if accessions.size == 1:
                self._identifier = accessions[0]

            return

        if self._name:

            loci = query.loc[:, 'Gene names  (primary )']
            loci_mask = []

            for value in loci:

                if self._name.lower() in value.lower() or value.lower() in self._name.lower():
                    loci_mask.append(True)

                else:
                    loci_mask.append(False)

            accessions = query.loc[loci_mask, 'Entry']

            if accessions.size == 1:
                self._identifier = accessions[0]

            return

    def fetch(self):

        if not self._identifier:

            if self._locus_tag and self._taxonomy:

                query = {'gene': self._locus_tag, 'taxonomy': self._taxonomy}

            elif self._locus_tag and not self._taxonomy:

                query = {'gene': self._locus_tag, 'taxonomy': self._taxonomy}

            elif self._name and self._taxonomy:

                query = {'gene': self._name, 'taxonomy': self._taxonomy}

            else:
                query = {}

            if query:
                uniprot_query = query_uniprot(query=query)

                # it sets up the uniprot accession
                self.parse_uniprot_query(uniprot_query)

        if not self._identifier:
            return

        self.record = fetch_uniprot_record(self._identifier, 'xml')
