import os
import shutil

import pandas as pd
import whoosh.index as w_index
from whoosh import searching
from whoosh.fields import Schema, TEXT
from whoosh.qparser import QueryParser

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

    cdw = os.getcwd()
    db_file = os.path.join(cdw, f'{db}.txt')
    url = f'{KEGGAPI.list_api}/{db}'

    if cache and os.path.exists(db_file):
        df_list = pd.read_csv(db_file, sep='\t', header=None)

    else:
        response = request(url)
        df_list = read_response(response, sep='\t', header=None)
        df_list.to_csv(f'{db}.txt', sep='\t', index=False)

    df_list.columns = ('identifier', 'name')

    return df_list


def indexing_kegg_list(df_kegg_list: pd.DataFrame, db: str, cache: bool = True) -> w_index.FileIndex:
    cdw = os.getcwd()
    index_dir = os.path.join(cdw, f'{db}_index')

    if cache and os.path.exists(index_dir):
        return w_index.open_dir(index_dir)

    elif not cache and os.path.exists(index_dir):
        shutil.rmtree(index_dir, ignore_errors=True)

    schema = Schema(identifier=TEXT(stored=True), name=TEXT(stored=True))

    index = w_index.create_in(index_dir, schema)
    writer = index.writer()

    for _, row in df_kegg_list.iterrows():
        writer.add_document(identifier=str(row['identifier']), name=str(row['name']))

    writer.commit()

    return index


def _results_to_set(results: searching.Results) -> set:

    results_set = set()
    for hit in results:

        hit_fields = hit.fields()
        kegg_db_id = hit_fields.get('identifier', '')
        kegg_id = kegg_db_id.split(':')

        if kegg_id:
            kegg_id = kegg_id[1]
            results_set.add(kegg_id)

    return results_set


def search_kegg_list(index: w_index.FileIndex, query: str, db: str) -> set:

    if db == 'pathway':
        forbidden_words = KEGGAPI.pathway_common_words

    else:
        forbidden_words = set()

    results = set()

    with index.searcher() as searcher:
        parser = QueryParser("name", index.schema)

        query_parser = parser.parse(query)
        query_results = searcher.search(query_parser)

        query_results = _results_to_set(query_results)
        results.update(query_results)

        if not results:

            for sub_query in query.split():

                if sub_query.lower() not in forbidden_words:

                    query_parser = parser.parse(sub_query)

                    sub_query_results = searcher.search(query_parser)

                    sub_query_results = _results_to_set(sub_query_results)
                    results.update(sub_query_results)

        return results
