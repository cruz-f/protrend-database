import pandas as pd

from protrend.transform.processors import remove_white_space, remove_regprecise_more, remove_multiple_white_space, \
    rstrip, lstrip, remove_pubmed, apply_processors
from protrend.transform.regprecise.settings import RegulatoryFamilySettings
from protrend.transform.transformer import Transformer


class RegulatoryFamilyTransformer(Transformer):

    def __init__(self, settings: RegulatoryFamilySettings = None):

        if not settings:
            settings = RegulatoryFamilySettings()

        super().__init__(settings)

    def read(self, **kwargs):
        return self._read_json_lines()

    @staticmethod
    def _transform_tf_family(df):
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

    @staticmethod
    def _transform_tf(df):

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

    @staticmethod
    def _transform_rna(df):
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

    def transform(self, **kwargs):

        tf_family = kwargs.get('tf_family', pd.DataFrame(columns=['name']))
        tf_family = self._transform_tf_family(tf_family)

        tf = kwargs.get('tf', pd.DataFrame(columns=['name']))
        tf = self._transform_tf(tf)

        rna = kwargs.get('rna', pd.DataFrame(columns=['name']))
        rna = self._transform_rna(rna)

        df = pd.merge(tf_family, tf, on='name', suffixes=('_tf_family', '_tf'))

        # concat description
        df['description'] = df['description_tf_family'].astype(str) + df['description_tf'].astype(str)

        # concat regulog
        df['regulog'] = df['regulog_tf_family'].astype(set) + df['regulog_tf'].astype(set)

        # concat pubmed
        df['pubmed'] = df['pubmed_tf_family'].astype(set) + df['pubmed_tf'].astype(set)

        df = df.drop(['description_tf_family', 'description_tf',
                      'regulog_tf_family', 'regulog_tf',
                      'pubmed_tf_family', 'pubmed_tf'], axis=1)

        df = pd.merge(df, rna, on='name', suffixes=('_tf', '_rna'))

        # concat description
        df['description'] = df['description_tf'].astype(str) + df['description_rna'].astype(str)

        # concat regulog
        df['regulog'] = df['regulog_tf'].astype(set) + df['regulog_rna'].astype(set)

        # concat pubmed
        df['pubmed'] = df['pubmed_tf'].astype(set) + df['pubmed_rna'].astype(set)

        df = df.drop(['description_tf', 'description_rna',
                      'regulog_tf', 'regulog_rna',
                      'pubmed_tf', 'pubmed_rna'],
                     axis=1)

        df = df.rename(columns={'url': 'url_rna'})

        source_dbs = []
        urls = []
        external_identifiers = []
        api_keys = []

        for _, row in df.iterrows():
            tf_family_id = row.get('tffamily_id', None)
            tf_id = row.get('collection_id', None)
            rna_id = row.get('riboswitch_id', None)

            if tf_family_id:
                db = 'regprecise'
                url = row.get('url_tf_family')
                external_identifier = tf_family_id
                api_key = 'tffamily_id'

            elif tf_id:
                db = 'regprecise'
                url = row.get('url_tf')
                external_identifier = tf_id
                api_key = 'collection_id'

            elif rna_id:
                db = 'regprecise'
                url = row.get('url_rna')
                external_identifier = rna_id
                api_key = 'riboswitch_id'

            else:
                db = None
                url = None
                external_identifier = None
                api_key = None

            source_dbs.append(db)
            urls.append(url)
            external_identifiers.append(external_identifier)
            api_keys.append(api_key)

        df['source_db'] = source_dbs
        df['url'] = urls
        df['external_identifier'] = external_identifiers
        df['api_key'] = api_keys

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df
