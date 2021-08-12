import os
from typing import Union, Tuple, Any

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from protrend.bioapis.settings import ENTREZ_E_MAIL, ENTREZ_API_KEY, ENTREZ_TOOL
from protrend.bioapis.utils import sleep, slugify
from protrend.io.json import read_json, write_json
from protrend.utils.settings import DATA_LAKE_BIOAPI_PATH

Entrez.email = ENTREZ_E_MAIL
Entrez.api_key = ENTREZ_API_KEY
Entrez.tool = ENTREZ_TOOL

ENTREZ_PATH = DATA_LAKE_BIOAPI_PATH.joinpath('entrez')
ENTREZ_SEARCH_PATH = DATA_LAKE_BIOAPI_PATH.joinpath('entrez', 'search')
ENTREZ_SUMMARY_PATH = DATA_LAKE_BIOAPI_PATH.joinpath('entrez', 'summary')
ENTREZ_FETCH_PATH = DATA_LAKE_BIOAPI_PATH.joinpath('entrez', 'fetch')

if not os.path.exists(ENTREZ_SEARCH_PATH):
    os.makedirs(ENTREZ_SEARCH_PATH)

if not os.path.exists(ENTREZ_SUMMARY_PATH):
    os.makedirs(ENTREZ_SUMMARY_PATH)

if not os.path.exists(ENTREZ_FETCH_PATH):
    os.makedirs(ENTREZ_FETCH_PATH)


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


def slugify_entrez(*args: Any) -> str:
    term = '/'.join(str(arg) for arg in args)
    term = slugify(term)
    return term


@sleep()
def entrez_search(db: str, term: str, retmax: int = 5000, cache: bool = True) -> Tuple[dict, bool]:
    cached_result = False
    slugified_term = slugify_entrez(db, term, retmax)
    file_path = os.path.join(ENTREZ_SEARCH_PATH, f'{slugified_term}.json')

    if os.path.exists(file_path) and cache:
        cached_result = True

        try:
            record = read_json(file_path)
        except:
            record = {}

    elif not os.path.exists(file_path) and cache:
        record = entrez_record(Entrez.esearch, db=db, term=term, retmax=retmax)

        try:
            write_json(file_path, record)
        except:
            record = {}

    else:
        record = entrez_record(Entrez.esearch, db=db, term=term, retmax=retmax)

    return record, cached_result


@sleep()
def entrez_summary(db: str, identifier: str, cache: bool = True) -> Tuple[dict, bool]:
    cached_result = False
    slugified_term = slugify_entrez(db, identifier)
    file_path = os.path.join(ENTREZ_SUMMARY_PATH, f'{slugified_term}.json')

    if os.path.exists(file_path) and cache:
        cached_result = True

        try:
            record = read_json(file_path)
        except:
            record = {}

    elif not os.path.exists(file_path) and cache:

        record = {}

        records = entrez_record(Entrez.esummary, db=db, id=identifier)

        if isinstance(records, list):
            record = records[0]

        elif isinstance(records, dict):
            d_sum_set = records.get('DocumentSummarySet', {})
            sum_set = d_sum_set.get('DocumentSummary', [{}])
            record = sum_set[0]

        try:
            write_json(file_path, record)
        except:
            record = {}

    else:
        record = {}

        records = entrez_record(Entrez.esummary, db=db, id=identifier)

        if isinstance(records, list):
            record = records[0]

        elif isinstance(records, dict):
            d_sum_set = records.get('DocumentSummarySet', {})
            sum_set = d_sum_set.get('DocumentSummary', [{}])
            record = sum_set[0]

    return record, cached_result


@sleep()
def entrez_fetch(db: str,
                 identifier: str,
                 rettype: str = "gb",
                 retmode: str = "text",
                 cache: bool = True) -> Tuple[SeqRecord, bool]:
    cached_result = False
    slugified_term = slugify_entrez(db, identifier, rettype, retmode)
    file_path = os.path.join(ENTREZ_FETCH_PATH, f'{slugified_term}.{rettype}')

    if rettype == 'gb':
        ret_type = 'genbank'

    else:
        ret_type = 'fasta'

    if os.path.exists(file_path) and cache:
        cached_result = True

        try:
            record = SeqIO.read(file_path, ret_type)

        except:
            record = SeqRecord(Seq(""))

    elif not os.path.exists(file_path) and cache:
        try:
            with Entrez.efetch(db=db, id=identifier, rettype=rettype, retmode=retmode) as handle:
                with open(file_path, "w") as file:
                    file.write(handle.read())
            record = SeqIO.read(file_path, ret_type)

        except:
            record = SeqRecord(Seq(""))

    else:
        try:
            handle = Entrez.efetch(db=db, id=identifier, rettype=rettype, retmode=retmode)
            record = SeqIO.read(handle, ret_type)

        except:
            record = SeqRecord(Seq(""))

    return record, cached_result
