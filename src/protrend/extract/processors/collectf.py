import re
from typing import Union


class CollecTFProcessors:

    pubmed_pattern = r'\[PMID:[0-9]*\]'

    @staticmethod
    def process_tax_onclick(onclick: str) -> Union[int, None]:
        if "getMotifReportsByTaxonomy" in onclick:
            tax_id = onclick.replace('getMotifReportsByTaxonomy(', '').replace('); return false;', '')
            return int(tax_id)

    @staticmethod
    def process_site_identifier(site_prop: str) -> str:

        if site_prop.lower() == 'act':
            return 'act'

        elif site_prop.lower() == 'rep':
            return 'rep'

        elif site_prop.lower() == 'n/a':
            return 'na'

        return str(site_prop)

    @staticmethod
    def process_mode(mode: str) -> str:

        if mode.lower() == 'act':
            return 'act'

        elif mode.lower() == 'rep':
            return 'rep'

        elif mode.lower() == 'n/a':
            return 'na'

        return 'na'

    @staticmethod
    def process_evidence_identifier(evidence_id: str) -> str:
        return re.sub(CollecTFProcessors.pubmed_pattern, '', evidence_id).strip()

    @staticmethod
    def process_pubmed(pubmed: str) -> Union[str, None]:
        match = re.search(CollecTFProcessors.pubmed_pattern, pubmed)

        if match:
            return match.group().replace('[PMID:', '').replace(']', '')

    @staticmethod
    def process_tf_pubmed(pubmed: str) -> str:
        return pubmed.replace('[', '').replace(']', '')