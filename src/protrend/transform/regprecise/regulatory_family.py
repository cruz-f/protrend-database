import pandas as pd

from protrend.io.csv import read_csv
from protrend.io.json import read_json_lines
from protrend.model.model import Source, Publication
from protrend.transform.processors import remove_white_space, remove_regprecise_more, remove_multiple_white_space, \
    rstrip, lstrip, remove_pubmed, apply_processors
from protrend.transform.regprecise.settings import RegulatoryFamilySettings
from protrend.transform.transformer import Transformer


class RegulatoryFamilyTransformer(Transformer):

    def __init__(self, settings: RegulatoryFamilySettings = None):

        if not settings:
            settings = RegulatoryFamilySettings()

        super().__init__(settings)

    def _transform_tf_family(self):
        file_path = self._transform_stack.get('tf_family')

        if not file_path:
            return pd.DataFrame(columns=['pubmed'])

        df = read_json_lines(file_path)

        df = df.drop_duplicates(subset=['name'])

        apply_processors(remove_white_space,
                         df=df,
                         col='name')

        apply_processors(remove_regprecise_more,
                         remove_pubmed,
                         remove_multiple_white_space,
                         rstrip,
                         lstrip,
                         df=df,
                         col='description')

        return df

    def _transform_tf(self):

        file_path = self._transform_stack.get('tf')

        if not file_path:
            return pd.DataFrame(columns=['pubmed'])

        df = read_json_lines(file_path)

        df = df.drop_duplicates(subset=['name'])

        apply_processors(remove_white_space,
                         df=df,
                         col='name')

        apply_processors(remove_regprecise_more,
                         remove_pubmed,
                         remove_multiple_white_space,
                         rstrip,
                         lstrip,
                         df=df,
                         col='description')

        return df

    def _transform_rna(self):
        file_path = self._transform_stack.get('rna')

        if not file_path:
            return pd.DataFrame(columns=['pubmed'])

        df = read_json_lines(file_path)

        df = df.drop_duplicates(subset=['name'])

        apply_processors(remove_white_space,
                         df=df,
                         col='name')

        apply_processors(remove_regprecise_more,
                         remove_pubmed,
                         remove_multiple_white_space,
                         rstrip,
                         lstrip,
                         df=df,
                         col='description')

        df = df.rename(columns={'url': 'url_rna'})

        # set mechanism
        size, _ = df.shape
        df['mechanism'] = ['small RNA (sRNA)'] * size

        return df

    @staticmethod
    def _transform_tfs(tf_family: pd.DataFrame, tf: pd.DataFrame) -> pd.DataFrame:

        df = pd.merge(tf_family, tf, on='name', suffixes=('_tf_family', '_tf'))

        # set mechanism
        size, _ = df.shape
        df['mechanism'] = ['transcription factor'] * size

        # concat description
        df['description'] = df['description_tf_family'].astype(str) + df['description_tf'].astype(str)

        # concat regulog
        df['regulog'] = df['regulog_tf_family'].astype(set) + df['regulog_tf'].astype(set)

        # concat pubmed
        df['pubmed'] = df['pubmed_tf_family'].astype(set) + df['pubmed_tf'].astype(set)

        df = df.drop(['description_tf_family', 'description_tf',
                      'regulog_tf_family', 'regulog_tf',
                      'pubmed_tf_family', 'pubmed_tf'], axis=1)

        return df

    def transform(self):

        # -------------------- TFs -------------------
        tf_family = self._transform_tf_family()
        tf = self._transform_tf()
        tfs = self._transform_tfs(tf_family=tf_family, tf=tf)

        rna = self._transform_rna()

        df = pd.concat([tfs, rna], axis=0)

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def _connect_to_source(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to_source')

        if not from_path:
            return pd.DataFrame()

        if not to_path:
            return pd.DataFrame()

        to_df = read_csv(to_path)
        to_df = to_df.query('name == regprecise')
        regprecise_id = to_df['protrend_id'].iloc[0]

        from_df = read_csv(from_path)

        tf_chain = [('tffamily_id', 'url_tf_family'),
                    ('collection_id', 'url_tf'),
                    ('riboswitch_id', 'url_rna')]

        dfs = []
        for key, url in tf_chain:
            key_mask = from_df[key].notnull()
            key_df = from_df[key_mask]

            from_identifiers = key_df['protrend_id'].tolist()
            size = len(from_identifiers)

            to_identifiers = [regprecise_id] * size

            kwargs = dict(url=key_df[url].tolist(),
                          external_identifier=key_df[key].tolist(),
                          key=[key] * size)

            df = self.make_connection(size=size,
                                      from_node=self.node,
                                      to_node=Source,
                                      from_identifiers=from_identifiers,
                                      to_identifiers=to_identifiers,
                                      kwargs=kwargs)

            dfs.append(df)

        return pd.concat(dfs, axis=0)

    def _connect_to_publication(self) -> pd.DataFrame:

        from_path = self._connect_stack.get('from')
        to_path = self._connect_stack.get('to_publication')

        if not from_path:
            return pd.DataFrame()

        if not to_path:
            return pd.DataFrame()

        from_df = read_csv(from_path)
        to_df = read_csv(to_path)

        from_identifiers = []
        to_identifiers = []

        for i, pubmed_row in enumerate(from_df['pubmed']):

            if not pubmed_row:
                continue

            for pmid in pubmed_row:

                from_id = from_df['protrend_id'].iloc[i]
                from_identifiers.append(from_id)

                mask = to_df['pmid'].values == pmid.replace(' ', '')
                to_id = to_df.loc[mask, 'protrend_id'].iloc[0]
                to_identifiers.append(to_id)

        size = len(from_identifiers)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Publication,
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers)

    def connect(self):

        source_connection = self._connect_to_source()

        df_name = f'connected_{self.node.node_name()}_{Source.node_name()}'
        self.stack_csv(df_name, source_connection)

        publication_connection = self._connect_to_publication()

        df_name = f'connected_{self.node.node_name()}_{Publication.node_name()}'
        self.stack_csv(df_name, publication_connection)
