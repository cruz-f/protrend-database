import io
from typing import Dict, Tuple, Union

import pandas as pd
from Bio import SeqIO
from Bio.SeqIO import SeqRecord

from protrend.utils.api_requests import request, read_response


class UniProtAPI:
    api = "https://www.uniprot.org/uniprot"
    mapping = "https://www.uniprot.org/uploadlists"
    query_fields = ('accession', 'ec', 'gene', 'gene_exact', 'id', 'organism', 'taxonomy')
    formats = ('html', 'tab', 'xls', 'fasta', 'gff', 'txt', 'xml', 'rdf', 'list', 'rss')
    column_names = ('id', 'entry name', 'genes', 'genes(PREFERRED)', 'genes(ALTERNATIVE)', 'genes(OLN)', 'organism',
                    'organism-id', 'protein names')
    custom_columns = ('id', 'entry name', 'genes', 'genes(PREFERRED)', 'organism', 'organism-id')
    custom_df_columns = ('Entry', 'Entry name', 'Gene names', 'Gene names  (primary )', 'Organism', 'Organism ID')
    mapping_terms = ('P_GI', 'P_ENTREZGENEID', 'REFSEQ_NT_ID', 'P_REFSEQ_AC', 'EMBL', 'EMBL_ID',
                     'ACC', 'ACC+ID', 'SWISSPROT', 'ID')


def fetch_uniprot_record(uniprot_accession: str, format_: str = 'xml') -> Union[None, SeqRecord, str]:
    if format_ not in UniProtAPI.formats:
        format_ = 'xml'

    url = f'{UniProtAPI.api}/{uniprot_accession}.{format_}'

    response = request(url)

    if not response.text:
        return

    if format_ == 'xml':
        handle = io.StringIO(response.text)
        return SeqIO.read(handle, 'uniprot-xml')

    return response.text


def query_uniprot(query: Dict[str, str],
                  format_: str = 'tab',
                  columns: Tuple[str] = None,
                  limit: str = 5,
                  sort: bool = True,
                  output: str = 'dataframe') -> Union[None, Dict, pd.DataFrame]:
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

    if format_ not in UniProtAPI.formats:
        format_ = 'tab'

    if not columns:
        columns = UniProtAPI.custom_columns

    else:
        columns = [col for col in columns if col in UniProtAPI.column_names]

        if not columns:
            columns = UniProtAPI.custom_columns

    columns_str = ','.join(columns)

    if sort:
        url = f'{UniProtAPI.api}/?query={query_str}&format={format_}&columns={columns_str}&limit={limit}&sort=score'

    else:
        url = f'{UniProtAPI.api}/?query={query_str}&format={format_}&columns={columns_str}&limit={limit}'

    response = request(url)

    if format_ == 'tab':
        sep = '\t'

        if output == 'dataframe':
            return read_response(response, sep=sep)

        elif output == 'dict':
            return read_response(response, sep=sep).to_dict()

    return response


def map_uniprot_identifiers(identifiers: Tuple[str],
                            from_: str,
                            to: str,
                            format_: str = 'tab',
                            output: str = 'dataframe') -> Union[None, Dict, pd.DataFrame]:

    if from_ not in UniProtAPI.mapping_terms:
        return

    if to not in UniProtAPI.mapping_terms:
        return

    query = ' '.join(identifiers)

    params = {'from': from_,
              'to': to,
              'format': format_,
              'query': query}

    url = f'{UniProtAPI.mapping}/'

    response = request(url, params=params)

    if format_ == 'tab':
        sep = '\t'

        if output == 'dataframe':
            return read_response(response, sep=sep)

        elif output == 'dict':
            return read_response(response, sep=sep).to_dict()

    return response
