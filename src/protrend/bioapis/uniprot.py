import io
import os
from typing import Dict, Tuple, Union, List

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord

from protrend.bioapis.utils import sleep, slugify
from protrend.utils.request import request, read_response
from protrend.utils.settings import DATA_LAKE_BIOAPI_PATH

UNIPROT_PATH = DATA_LAKE_BIOAPI_PATH.joinpath('uniprot')
UNIPROT_RECORDS_PATH = DATA_LAKE_BIOAPI_PATH.joinpath('uniprot', 'records')
UNIPROT_QUERY_PATH = DATA_LAKE_BIOAPI_PATH.joinpath('uniprot', 'query')

if not os.path.exists(UNIPROT_RECORDS_PATH):
    os.makedirs(UNIPROT_RECORDS_PATH)

if not os.path.exists(UNIPROT_QUERY_PATH):
    os.makedirs(UNIPROT_QUERY_PATH)


class UniProtAPI:
    record = "https://www.uniprot.org/uniprot"
    record_format = 'xml'
    query_fields = ('accession', 'ec', 'gene', 'gene_exact', 'id', 'organism', 'taxonomy')
    query_columns = ('id', 'entry name', 'genes', 'genes(PREFERRED)', 'organism', 'organism-id')
    query_format = 'tab'
    df_query_columns = ('Entry', 'Entry name', 'Gene names', 'Gene names  (primary )', 'Organism', 'Organism ID')

    mapping = "https://www.uniprot.org/uploadlists"
    mapping_format = 'tab'
    mapping_terms = ('P_GI', 'P_ENTREZGENEID', 'REFSEQ_NT_ID', 'P_REFSEQ_AC', 'EMBL', 'EMBL_ID',
                     'ACC', 'ACC+ID', 'SWISSPROT', 'ID')
    df_mapping_columns = ('From', 'To')


def slugify_uniprot(url: str) -> str:
    columns_str = ','.join(UniProtAPI.query_columns)
    url = slugify(url).replace('httpswwwuniprotorguniprotquery', '')
    url = url.replace(f'&format={UniProtAPI.query_format}&columns={columns_str}', '')
    return url


@sleep(0.5)
def fetch_uniprot_record(uniprot_accession: str, cache: bool = True) -> Tuple[SeqRecord, bool]:
    cached_result = False
    url = f'{UniProtAPI.record}/{uniprot_accession}.{UniProtAPI.record_format}'
    file_path = os.path.join(UNIPROT_RECORDS_PATH, f'{uniprot_accession}.xml')

    if os.path.exists(file_path) and cache:
        handle = open(file_path, 'r')
        cached_result = True

    elif not os.path.exists(file_path) and cache:
        response = request(url)

        with open(file_path, 'w') as file:
            file.write(response.text)

        handle = open(file_path, 'r')

    else:
        response = request(url)
        handle = io.StringIO(response.text)

    try:
        return SeqIO.read(handle, 'uniprot-xml'), cached_result

    except:
        return SeqRecord(Seq("")), cached_result


@sleep(0.5)
def query_uniprot(query: Dict[str, str],
                  limit: int = 5,
                  sort: bool = True,
                  cache: bool = True,
                  output: str = 'dataframe') -> Tuple[Union[Dict, pd.DataFrame], bool]:
    cached_result = False
    query_str = ''

    for key, val in query.items():

        if key in UniProtAPI.query_fields:
            query_str += f'{key}:{val}+AND+'

        else:
            query_str += f'{val}+AND+'

    if not query_str:
        cached_result = True
        df = pd.DataFrame(columns=UniProtAPI.df_query_columns)

        if output == 'dataframe':
            return df, cached_result

        return df.to_dict(), cached_result

    # remove +AND+
    query_str = query_str[:-5]

    columns_str = ','.join(UniProtAPI.query_columns)

    if sort:
        url = f'{UniProtAPI.record}' \
              f'/?query={query_str}&format={UniProtAPI.query_format}&columns={columns_str}&limit={limit}&sort=score'

    else:
        url = f'{UniProtAPI.record}' \
              f'/?query={query_str}&format={UniProtAPI.query_format}&columns={columns_str}&limit={limit}'

    slugified_url = slugify_uniprot(url)

    file_path = os.path.join(UNIPROT_RECORDS_PATH, f'{slugified_url}.tsv')

    if os.path.exists(file_path) and cache:
        cached_result = True
        df = pd.read_csv(file_path, sep='\t')

    elif not os.path.exists(file_path) and cache:
        response = request(url)
        df = read_response(response, sep='\t')

        df.to_csv(file_path, sep='\t', index=False)

    else:
        response = request(url)
        df = read_response(response, sep='\t')

    if df.empty:
        df = pd.DataFrame(columns=UniProtAPI.df_query_columns)

    if output == 'dataframe':
        return df, cached_result

    return df.to_dict(), cached_result


@sleep(0.5)
def map_uniprot_identifiers(identifiers: Union[List[str], Tuple[str]],
                            from_: str,
                            to: str,
                            output: str = 'dataframe') -> Union[Dict, pd.DataFrame]:
    if from_ not in UniProtAPI.mapping_terms:
        raise ValueError(f'Invalid from {from_}')

    if to not in UniProtAPI.mapping_terms:
        raise ValueError(f'Invalid from {to}')

    query = ' '.join(identifiers)

    params = {'from': from_,
              'to': to,
              'format': UniProtAPI.mapping_format,
              'query': query}

    url = f'{UniProtAPI.mapping}/'

    response = request(url, params=params)

    df = read_response(response, sep='\t')

    if df.empty:
        df = pd.DataFrame(columns=UniProtAPI.df_mapping_columns)

    if output == 'dataframe':
        return df

    return df.to_dict()
