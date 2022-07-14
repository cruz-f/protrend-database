from typing import List

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.entrez import entrez_summary, entrez_search, entrez_fetch
from protrend.bioapis.uniprot import fetch_uniprot_record, query_uniprot
from protrend.utils.processors import to_int_str, to_str, apply_processors, lower_case, split_str
from protrend.utils.miscellaneous import is_null


class NCBIProtein(BioAPI):

    def __init__(self,
                 identifier: str = '',
                 refseq_accession: str = '',
                 genbank_accession: str = '',
                 taxonomy: str = '',
                 locus_tag: str = '',
                 name: str = ''):

        super().__init__(identifier)

        if is_null(refseq_accession):
            refseq_accession = ''

        refseq_accession = str(refseq_accession)

        if is_null(genbank_accession):
            genbank_accession = ''

        genbank_accession = str(genbank_accession)

        if is_null(taxonomy):
            taxonomy = ''

        taxonomy = to_int_str(taxonomy)

        if is_null(locus_tag):
            locus_tag = ''

        locus_tag = str(locus_tag)

        if is_null(name):
            name = ''

        name = str(name)

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
    def taxonomy(self) -> str:
        return str(self.record.get('TaxId', self._taxonomy))

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
    def seq_record(self) -> SeqRecord:
        return self._seq_record

    @seq_record.setter
    def seq_record(self, value):
        if value:
            self._seq_record = value

    @property
    def synonyms(self) -> List[str]:

        synonyms = []

        if self._locus_tag:
            synonyms.append(self._locus_tag)

        if self._name:
            synonyms.append(self._name)

        if self.locus_tag:
            synonyms.append(self.locus_tag)

        if self.seq_record:
            for feature in self.seq_record.features:

                if feature.type == 'CDS':
                    old_locus_tag = feature.qualifiers.get('old_locus_tag')

                    if old_locus_tag:
                        synonyms.append(old_locus_tag[0])

                        break

        return synonyms

    @property
    def sequence(self) -> str:
        return str(self.seq_record.seq)

    # noinspection DuplicatedCode
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

    def fetch(self, *args, **kwargs):

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
            record = entrez_summary(db='protein', identifier=identifier)
            self.record = record

            record = entrez_fetch(db='protein', identifier=identifier)
            self.seq_record = record


class UniProtProtein(BioAPI):

    def __init__(self,
                 identifier: str = '',
                 taxonomy: str = '',
                 locus_tag: str = '',
                 name: str = ''):

        super().__init__(identifier)

        if is_null(taxonomy):
            taxonomy = ''

        taxonomy = to_int_str(taxonomy)

        # E. coli taxonomy identifier for uniprot is 83333 rather than 511145
        if taxonomy == '511145':
            taxonomy = '83333'

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
        return getattr(self.record, 'id', self._identifier)

    @property
    def taxonomy(self) -> str:
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

        if hasattr(self.record, 'annotations'):

            if 'gene_name_primary' in self.record.annotations:
                return self.record.annotations['gene_name_primary']

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

        if self._locus_tag:
            synonyms.append(self._locus_tag)

        if self._name:
            synonyms.append(self._name)

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

    def _filter_by_taxonomy(self, query: pd.DataFrame) -> pd.DataFrame:

        if self._taxonomy:
            tax_mask = query['Organism (ID)'].str.contains(self._taxonomy, na=False, regex=False)

            if tax_mask.any():
                query = query[tax_mask]

        return query

    def _filter_by_locus_tag(self, query: pd.DataFrame) -> pd.DataFrame:

        def loci_filter(row):

            if is_null(row):
                return False

            return self._locus_tag.lower() in row

        if self._locus_tag:
            loci_mask = query['Gene Names'].map(loci_filter)

            if loci_mask.any():
                query = query[loci_mask]

        return query

    def _filter_by_name(self, query: pd.DataFrame) -> pd.DataFrame:

        def name_filter(row):

            if is_null(row):
                return False

            return self._name.lower() in row

        if self._name:
            name_mask = query['Gene Names (primary)'].map(name_filter)

            if name_mask.any():
                query = query[name_mask]

        return query

    def parse_uniprot_query(self, query: pd.DataFrame):

        processors = {'Organism (ID)': to_int_str,
                      'Gene Names': [to_str, lower_case, split_str],
                      'Gene Names (primary)': [to_str, lower_case, split_str]}
        query = apply_processors(query, **processors)

        query = self._filter_by_taxonomy(query=query)
        query = self._filter_by_locus_tag(query=query)
        query = self._filter_by_name(query=query)

        if query['Entry'].size == 1:
            return query['Entry'].iloc[0]

        return ''

    def build_query(self) -> dict:

        if self._locus_tag and self._taxonomy:

            return {'gene': self._locus_tag, 'taxonomy_id': self._taxonomy}

        elif self._locus_tag and not self._taxonomy:

            return {'gene': self._locus_tag}

        elif self._name and self._taxonomy:

            return {'gene': self._name, 'taxonomy_id': self._taxonomy}

        return {}

    def fetch(self, *args, **kwargs):

        if self._identifier:
            identifier = self._identifier

        else:

            query = self.build_query()
            uniprot_query = query_uniprot(query=query)
            identifier = self.parse_uniprot_query(uniprot_query)

        if identifier:
            self._identifier = identifier

        else:
            return

        record = fetch_uniprot_record(self._identifier)
        if record.id == '<unknown id>':
            record = {}
        self.record = record
