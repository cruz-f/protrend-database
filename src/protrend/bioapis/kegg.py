import os
import shutil
from typing import Tuple, Set

import pandas as pd
import whoosh.index as w_index
from whoosh import searching
from whoosh.fields import Schema, TEXT
from whoosh.qparser import QueryParser

from protrend.bioapis.utils import BIO_APIS_DIR
from protrend.utils.api_requests import request, read_response


class KEGGAPI:
    list_api = 'http://rest.kegg.jp/list'

    dbs = ('pathway', 'compound')

    pathway_common_words = {'biosynthesis',
                            'pathway',
                            'signaling',
                            'metabolism',
                            'degradation',
                            'receptor',
                            'acid',
                            'cancer',
                            'agents',
                            'cell',
                            'antagonists',
                            'infection',
                            'secretion',
                            'agonists',
                            'disease',
                            'resistance',
                            'drugs',
                            'inhibitors',
                            'hormone',
                            'type',
                            'transport',
                            'cycle',
                            'protein',
                            'yeast',
                            'pathways',
                            'plant',
                            'carbon',
                            'system'}


def fetch_kegg_list(db: str, cache: bool = True) -> pd.DataFrame:

    if db not in KEGGAPI.dbs:
        raise ValueError(f'Invalid KEGG database: {db}')

    db_file = os.path.join(BIO_APIS_DIR, db, f'{db}.txt')
    url = f'{KEGGAPI.list_api}/{db}'

    if cache and os.path.exists(db_file):
        df_list = pd.read_csv(db_file, sep='\t', header=None)

    else:
        response = request(url)
        df_list = read_response(response, sep='\t', header=None)
        df_list.to_csv(db_file, sep='\t', index=False)

    df_list.columns = ('identifier', 'name')

    return df_list


def indexing_kegg_list(db: str, df_kegg_list: pd.DataFrame = None, cache: bool = True) -> w_index.FileIndex:
    index_dir = os.path.join(BIO_APIS_DIR, f'{db}_index')

    if cache and os.path.exists(index_dir):
        return w_index.open_dir(index_dir)

    elif not cache and os.path.exists(index_dir):
        shutil.rmtree(index_dir, ignore_errors=True)

    if df_kegg_list is None:
        raise ValueError(f'Invalid input: Index was not found in the default directory {index_dir} '
                         f'and kegg list is empty')

    schema = Schema(identifier=TEXT(stored=True), name=TEXT(stored=True))

    index = w_index.create_in(index_dir, schema)
    writer = index.writer()

    for _, row in df_kegg_list.iterrows():
        writer.add_document(identifier=str(row['identifier']), name=str(row['name']))

    writer.commit()

    return index


def _identifiers_to_set(results: searching.Results) -> set:
    results_set = set()
    for hit in results:

        hit_fields = hit.fields()
        kegg_db_id = hit_fields.get('identifier', '')
        kegg_id = kegg_db_id.split(':')

        if kegg_id:
            kegg_id = kegg_id[1]
            results_set.add(kegg_id)

    return results_set


def _names_to_set(results):

    results_set = set()
    for hit in results:

        hit_fields = hit.fields()
        kegg_db_names = hit_fields.get('name', '')
        kegg_names = kegg_db_names.split(';')

        if len(kegg_names) == 1:
            results_set.add(kegg_names[0])

        elif len(kegg_names) > 1:
            results_set.update(kegg_names)

    return results_set


def _results_to_set(results: searching.Results, field: str) -> set:

    if field == 'identifier':
        return _identifiers_to_set(results)

    return _names_to_set(results)


def search_kegg_list(index: w_index.FileIndex, query: str, db: str) -> Tuple[Set[str], Set[str]]:

    if db == 'pathway':
        forbidden_words = KEGGAPI.pathway_common_words

    else:
        forbidden_words = set()

    identifiers_results = set()
    names_results = set()

    with index.searcher() as searcher:
        parser = QueryParser("name", index.schema)

        query_parser = parser.parse(query)
        query_results = searcher.search(query_parser)

        identifiers_results.update(_results_to_set(query_results, 'identifier'))
        names_results.update(_results_to_set(query_results, 'name'))

        if not identifiers_results:

            for sub_query in query.split():

                if sub_query.lower() not in forbidden_words:

                    query_parser = parser.parse(sub_query)

                    sub_query_results = searcher.search(query_parser)

                    identifiers_results.update(_results_to_set(sub_query_results, 'identifier'))
                    names_results.update(_results_to_set(query_results, 'name'))

        return identifiers_results, names_results
