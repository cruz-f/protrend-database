import io
import os
import time
from typing import Dict, Tuple, Union, List

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from diskcache import Cache

from protrend.bioapis._map_identifiers_utils import submit_id_mapping, check_id_mapping_status, \
    get_id_mapping_results_link, get_id_mapping_results_search
from protrend.log import ProtrendLogger
from protrend.utils import Settings, request, read_response


def _init_uniprot_record() -> Cache:
    if not os.path.exists(Settings.uniprot_record):
        os.makedirs(Settings.uniprot_record)

    cache = Cache(directory=Settings.uniprot_record)
    return cache


def _init_uniprot_query() -> Cache:
    if not os.path.exists(Settings.uniprot_query):
        os.makedirs(Settings.uniprot_query)

    cache = Cache(directory=Settings.uniprot_query)
    return cache


def _init_uniprot_mapping() -> Cache:
    if not os.path.exists(Settings.uniprot_mapping):
        os.makedirs(Settings.uniprot_mapping)

    cache = Cache(directory=Settings.uniprot_mapping)
    return cache


record_cache = _init_uniprot_record()
query_cache = _init_uniprot_query()
mapping_cache = _init_uniprot_mapping()


class UniProtAPI:
    record = "https://rest.uniprot.org/uniprotkb"
    record_format = 'xml'
    query_fields = ('accession', 'ec', 'gene', 'gene_exact', 'accession_id', 'organism_name', 'organism_id',
                    'taxonomy_id', 'taxonomy_name')
    query_columns = ('accession', 'id', 'gene_names', 'gene_primary', 'organism_name', 'organism_id')
    query_format = 'tsv'
    df_query_columns = ('Entry', 'Entry Name', 'Gene Names', 'Gene Names (primary)', 'Organism', 'Organism (ID)')


@record_cache.memoize()
def fetch_uniprot_record(uniprot_accession: str) -> SeqRecord:
    url = f'{UniProtAPI.record}/{uniprot_accession}.{UniProtAPI.record_format}'
    response = request(url)
    handle = io.StringIO(response.text)

    try:
        record = SeqIO.read(handle, 'uniprot-xml')

    except:
        record = SeqRecord(Seq(""))

    time.sleep(Settings.request_sleep)

    return record


@query_cache.memoize()
def query_uniprot(query: Dict[str, str],
                  size: int = 5,
                  output: str = 'dataframe') -> Union[Dict, pd.DataFrame]:
    query_str = ''

    for key, val in query.items():

        if key in UniProtAPI.query_fields:
            query_str += f'{key}:{val}+AND+'

        else:
            query_str += f'{val}+AND+'

    if not query_str:
        df = pd.DataFrame(columns=UniProtAPI.df_query_columns)

        if output == 'dataframe':
            return df

        return df.to_dict()

    # remove +AND+
    query_str = query_str[:-5]

    columns_str = ','.join(UniProtAPI.query_columns)

    url = f'{UniProtAPI.record}/search?query={query_str}&format={UniProtAPI.query_format}&fields={columns_str}&size={size}'

    response = request(url)
    df = read_response(response, sep='\t')

    if df.empty:
        df = pd.DataFrame(columns=UniProtAPI.df_query_columns)

    time.sleep(Settings.request_sleep)

    if output == 'dataframe':
        return df

    return df.to_dict()


@mapping_cache.memoize()
def map_uniprot_identifiers(identifiers: Union[List[str], Tuple[str]],
                            from_db: str,
                            to_db: str,
                            output: str = 'dataframe') -> Union[Dict, pd.DataFrame]:

    job_id = submit_id_mapping(from_db=from_db,
                               to_db=to_db,
                               ids=identifiers)

    try:
        if check_id_mapping_status(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)
        else:
            results = {}

    except Exception as e:
        ProtrendLogger.log.error(e)
        results = {}

    df = pd.DataFrame(results.get('results', []))

    if df.empty:
        df = pd.DataFrame(columns=['from', 'to'])

    time.sleep(Settings.request_sleep)

    if output == 'dataframe':
        return df

    return df.to_dict()
