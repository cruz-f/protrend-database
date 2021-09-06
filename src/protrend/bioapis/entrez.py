import os
import time
from typing import Union

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from diskcache import Cache, JSONDisk

from protrend.bioapis.settings import ENTREZ_E_MAIL, ENTREZ_API_KEY, ENTREZ_TOOL
from protrend.utils.settings import DATA_LAKE_BIOAPI_PATH, REQUEST_SLEEP

Entrez.email = ENTREZ_E_MAIL
Entrez.api_key = ENTREZ_API_KEY
Entrez.tool = ENTREZ_TOOL

ENTREZ_PATH = DATA_LAKE_BIOAPI_PATH.joinpath('entrez')


def _init_entrez_search() -> Cache:
    directory = ENTREZ_PATH.joinpath('search')

    if not os.path.exists(directory):
        os.makedirs(directory)

    cache = Cache(directory=directory, disk=JSONDisk)

    return cache


def _init_entrez_summary() -> Cache:
    directory = ENTREZ_PATH.joinpath('summary')

    if not os.path.exists(directory):
        os.makedirs(directory)

    cache = Cache(directory=directory, disk=JSONDisk)

    return cache


def _init_entrez_fetch() -> Cache:
    directory = ENTREZ_PATH.joinpath('fetch')

    if not os.path.exists(directory):
        os.makedirs(directory)

    cache = Cache(directory=directory)

    return cache


search_cache = _init_entrez_search()
summary_cache = _init_entrez_summary()
fetch_cache = _init_entrez_fetch()


def entrez_record(callback, **kwargs) -> Union[SeqRecord, dict]:
    handle = None

    try:
        handle = callback(**kwargs)
        return Entrez.read(handle)

    except:
        return {}

    finally:
        if hasattr(handle, 'close'):
            handle.close()


@search_cache.memoize()
def entrez_search(db: str, term: str, retmax: int = 5000) -> dict:
    record = entrez_record(Entrez.esearch, db=db, term=term, retmax=retmax)

    time.sleep(REQUEST_SLEEP)
    return record


@summary_cache.memoize()
def entrez_summary(db: str, identifier: str) -> dict:
    records = entrez_record(Entrez.esummary, db=db, id=identifier)

    if isinstance(records, list):
        record = records[0]

    elif isinstance(records, dict):
        d_sum_set = records.get('DocumentSummarySet', {})
        sum_set = d_sum_set.get('DocumentSummary', [{}])
        record = sum_set[0]

    else:
        record = {}

    if hasattr(record, 'attributes'):
        record['attributes'] = record.attributes

    time.sleep(REQUEST_SLEEP)

    return record


@fetch_cache.memoize()
def entrez_fetch(db: str,
                 identifier: str,
                 rettype: str = "gb",
                 retmode: str = "text") -> SeqRecord:

    if rettype == 'gb':
        ret_type = 'genbank'

    else:
        ret_type = 'fasta'

    try:
        with Entrez.efetch(db=db, id=identifier, rettype=rettype, retmode=retmode) as handle:
            record = SeqIO.read(handle, ret_type)

    except:
        record = SeqRecord(Seq(""))

    time.sleep(REQUEST_SLEEP)

    return record
