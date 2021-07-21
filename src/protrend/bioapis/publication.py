from protrend.bioapis.bioapi import BioAPI
from protrend.bioapis.entrez import entrez_summary


class Publication(BioAPI):

    @property
    def pmid(self) -> str:
        return self.identifier

    @property
    def doi(self) -> str:
        return self.record.get('DOI')

    @property
    def title(self) -> str:
        return self.record.get('Title')

    @property
    def author(self) -> str:
        return self.record.get('LastAuthor')

    @property
    def year(self) -> str:

        pub_date = self.record.get('PubDate')

        if pub_date:
            year = ''.join(char for char in pub_date if char.isdigit())

            if year:
                return year

    def fetch(self):

        if self.pmid:
            self.record = entrez_summary(db='pubmed', identifier=self.pmid)
