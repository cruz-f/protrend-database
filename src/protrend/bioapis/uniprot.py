import io
from typing import Dict, Tuple, Union

import pandas as pd
from Bio import SeqIO
from Bio.SeqIO import SeqRecord

from protrend.bioapis.utils import sleep
from protrend.utils.api_requests import request, read_response


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


@sleep(0.5)
def fetch_uniprot_record(uniprot_accession: str) -> Union[None, SeqRecord]:

    url = f'{UniProtAPI.record}/{uniprot_accession}.{UniProtAPI.record_format}'

    response = request(url)

    if not response.text:
        return

    handle = io.StringIO(response.text)
    return SeqIO.read(handle, 'uniprot-xml')


@sleep(0.5)
def query_uniprot(query: Dict[str, str],
                  limit: str = 5,
                  sort: bool = True,
                  output: str = 'dataframe') -> Union[Dict, pd.DataFrame]:

    query_str = ''

    for key, val in query.items():

        if key in UniProtAPI.query_fields:
            query_str += f'{key}:{val}+AND+'

        else:
            query_str += f'{val}+AND+'

    if not query_str:
        return

    # remove +AND+
    query_str = query_str[:-5]

    columns_str = ','.join(UniProtAPI.query_columns)

    if sort:
        url = f'{UniProtAPI.record}/?query={query_str}&format={UniProtAPI.query_format}&columns={columns_str}&limit={limit}&sort=score'

    else:
        url = f'{UniProtAPI.record}/?query={query_str}&format={UniProtAPI.query_format}&columns={columns_str}&limit={limit}'

    response = request(url)

    df = read_response(response, sep='\t')

    if df.empty:
        df = pd.DataFrame(columns=UniProtAPI.df_query_columns)

    if output == 'dataframe':

        return df

    return df.to_dict()


def map_uniprot_identifiers(identifiers: Tuple[str],
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
