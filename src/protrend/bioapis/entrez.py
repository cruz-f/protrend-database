from contextlib import contextmanager

from Bio import Entrez, SeqIO

from protrend.bioapis.exception import BioException
from protrend.bioapis.settings import ENTREZ_E_MAIL, ENTREZ_API_KEY, ENTREZ_TOOL
from protrend.bioapis.utils import sleep

Entrez.email = ENTREZ_E_MAIL
Entrez.api_key = ENTREZ_API_KEY
Entrez.tool = ENTREZ_TOOL


@contextmanager
def entrez_record(callback, **kwargs):

    handle = None

    try:
        handle = callback(**kwargs)

        if callback == Entrez.efetch:

            rettype = kwargs.get('rettype', 'genbank')

            if rettype == 'gb':
                rettype = 'genbank'

            else:
                rettype = 'fasta'

            yield SeqIO.read(handle, rettype)

        else:
            yield Entrez.read(handle)

    except BioException as exc:
        yield {}

    finally:
        if hasattr(handle, 'close'):
            handle.close()

@sleep()
def entrez_search(db, term, retmax=5000):
    with entrez_record(Entrez.esearch, db=db, term=term, retmax=retmax) as record:
        return record


@sleep()
def entrez_summary(db, identifier):
    with entrez_record(Entrez.esummary, db=db, id=identifier) as records:

        if records and isinstance(records, list):
            return records[0]

        elif records and isinstance(records, dict):
            d_sum_set = records.get('DocumentSummarySet', {})
            sum_set = d_sum_set.get('DocumentSummary', [{}])
            return sum_set[0]

        return {}


@sleep()
def entrez_fetch(db, identifier, rettype="gb", retmode="text"):
    with entrez_record(Entrez.efetch, db=db, id=identifier, rettype=rettype, retmode=retmode) as record:
        return record
