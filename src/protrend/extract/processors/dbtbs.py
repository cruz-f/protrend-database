from typing import Union


class DBTBSProcessors:

    @staticmethod
    def process_pubmed(href: str) -> Union[int, None]:
        # http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=PubMed&dopt=Abstract&list_uids=+12100558

        if "list_uids=+" in href:
            *_, identifier = href.split("list_uids=+")
            return int(identifier)

        return

    @staticmethod
    def process_nd(text: str) -> Union[str, None]:
        # ND

        if isinstance(text, str):

            if "nd" == text.lower():
                return

        return text
