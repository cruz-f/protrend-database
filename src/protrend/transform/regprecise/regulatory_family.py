import pandas as pd

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

    def read(self, **kwargs):
        files = self._read_json_lines()

        return super(RegulatoryFamilyTransformer, self).read(**files)

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

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df

    def _connect_to_source(self, df: pd.DataFrame, name: str) -> pd.DataFrame:
        source_mask = df[name].notnull()
        df = df.loc[source_mask, :]
        size, _ = df.shape

        from_identifiers = list(df['protrend_id'])
        to_identifiers = ['regprecise'] * size

        kwargs = dict(url=list(df['url']),
                      external_identifier=list(df[name]),
                      name=[name] * size)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Source,
                                    from_property=self.node.identifying_property,
                                    to_property='name',
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers,
                                    **kwargs)

    def _connect_to_publication(self, df: pd.DataFrame) -> pd.DataFrame:
        from_identifiers = []
        to_identifiers = []

        for _, row in df.iterrows():

            pubmed_row = row.get('pubmed', [])

            for identifier in pubmed_row:

                from_identifiers.append(df[self.node.identifying_property])
                to_identifiers.append(identifier)

        size = len(from_identifiers)

        return self.make_connection(size=size,
                                    from_node=self.node,
                                    to_node=Publication,
                                    from_property=self.node.identifying_property,
                                    to_property='pmid',
                                    from_identifiers=from_identifiers,
                                    to_identifiers=to_identifiers)

    def connect(self, df: pd.DataFrame):

        tf_family_connection = self._connect_to_source(df, 'tffamily_id')
        tf_connection = self._connect_to_source(df, 'collection_id')
        rna_connection = self._connect_to_source(df, 'riboswitch_id')

        source_connection = pd.concat([tf_family_connection, tf_connection, rna_connection], axis=0)

        df_name = f'connected_{self.node.node_name()}_{Source.node_name()}'
        self.stack_csv(df_name, source_connection)

        publication_connection = self._connect_to_publication(df)

        df_name = f'connected_{self.node.node_name()}_{Publication.node_name()}'
        self.stack_csv(df_name, publication_connection)
