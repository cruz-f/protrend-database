import pandas as pd

from protrend.io.utils import read_from_stack
from protrend.transform.connector import DefaultConnector
from protrend.transform.processors import (remove_white_space, remove_regprecise_more, remove_multiple_white_space,
                                           rstrip, lstrip, remove_pubmed, apply_processors)
from protrend.transform.regprecise.publication import PublicationTransformer
from protrend.transform.regprecise.settings import RegulatoryFamilySettings, RegulatoryFamilyToSource, \
    RegulatoryFamilyToPublication
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import DefaultTransformer


class RegulatoryFamilyTransformer(DefaultTransformer):
    default_settings = RegulatoryFamilySettings
    columns = {'protrend_id',
               'mechanism'
               'tffamily_id', 'riboswitch_id', 'collection_id', 'name', 'url_tf_family', 'url_tf', 'url_rna',
               'description', 'pubmed', 'rfam', 'regulog'}

    tf_family_columns = {'tffamily_id', 'name', 'url', 'description', 'pubmed', 'regulog'}
    tf_columns = {'collection_id', 'name', 'url', 'description', 'pubmed', 'regulog'}
    rna_columns = {'riboswitch_id', 'name', 'url', 'description', 'pubmed', 'rfam', 'regulog'}

    def _transform_tf_family(self, tf_family: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=tf_family, subset=['name'], perfect_match=True, preserve_nan=False)

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

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True, preserve_nan=False)

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

    def _transform_rna(self, rna: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=rna, subset=['rfam'], perfect_match=True, preserve_nan=False)

        apply_processors(remove_white_space,
                         df=df,
                         col='rfam')

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

        apply_processors(set, df=df, col='pubmed')
        apply_processors(set, df=df, col='regulog')

        df = df.rename(columns={'url': 'url_rna'})

        # set mechanism
        size, _ = df.shape
        df['mechanism'] = ['small RNA (sRNA)'] * size

        return df

    def _transform_tfs(self, tf_family: pd.DataFrame, tf: pd.DataFrame) -> pd.DataFrame:
        df = pd.merge(tf_family, tf, on='name', suffixes=('_tf_family', '_tf'))

        # set mechanism
        size, _ = df.shape
        df['mechanism'] = ['transcription factor'] * size

        # concat description
        df = self.merge_columns(df=df, column='description', left='description_tf_family',
                                right='description_tf', fill='')

        # concat pubmed
        apply_processors(set, df=df, col='pubmed_tf_family')
        apply_processors(set, df=df, col='pubmed_tf')
        df = self.merge_columns(df=df, column='pubmed', left='pubmed_tf_family',
                                right='regulog_tf', fill=set())

        # concat regulog
        apply_processors(set, df=df, col='regulog_tf_family')
        apply_processors(set, df=df, col='regulog_tf')
        df = self.merge_columns(df=df, column='regulog', left='regulog_tf_family',
                                right='regulog_tf', fill=set())

        return df

    def transform(self):
        # -------------------- TFs -------------------
        tf_family = read_from_stack(tl=self, file='tf_family', json=True, default_columns=self.tf_family_columns)
        tf_family = self._transform_tf_family(tf_family)
        tf = read_from_stack(tl=self, file='tf', json=True, default_columns=self.tf_columns)
        tf = self._transform_tf(tf)
        tfs = self._transform_tfs(tf_family=tf_family, tf=tf)

        # -------------------- RNAs -------------------
        rna = read_from_stack(tl=self, file='rna', json=True, default_columns=self.rna_columns)
        rna = self._transform_rna(rna)

        df = pd.concat([tfs, rna], axis=0)

        if df.empty:
            df = self.make_empty_frame()

        df_name = f'transformed_{self.node.node_name()}'
        self.stack_csv(df_name, df)

        return df


class RegulatoryFamilyToSourceConnector(DefaultConnector):
    default_settings = RegulatoryFamilyToSource

    def connect(self):
        regulatory_family = read_from_stack(tl=self, file='regulatory_family', json=False,
                                            default_columns=RegulatoryFamilyTransformer.columns)
        source = read_from_stack(tl=self, file='source', json=False, default_columns=SourceTransformer.columns)

        protrend_id = source['protrend_id'].iloc[0]

        tf_chain = [('tffamily_id', 'url_tf_family'),
                    ('collection_id', 'url_tf'),
                    ('riboswitch_id', 'url_rna')]

        dfs = []
        for key, url in tf_chain:
            key_mask = regulatory_family[key].notnull()
            key_df = regulatory_family[key_mask]

            from_identifiers = key_df['protrend_id'].tolist()
            size = len(from_identifiers)

            to_identifiers = [protrend_id] * size

            kwargs = dict(url=key_df[url].tolist(),
                          external_identifier=key_df[key].tolist(),
                          key=[key] * size)

            df = self.make_connection(size=size,
                                      from_identifiers=from_identifiers,
                                      to_identifiers=to_identifiers,
                                      kwargs=kwargs)

            dfs.append(df)

        df = pd.concat(dfs, axis=0)

        self.stack_csv(df)


class RegulatoryFamilyToPublicationConnector(DefaultConnector):
    default_settings = RegulatoryFamilyToPublication

    def connect(self):
        regulatory_family = read_from_stack(tl=self, file='regulatory_family', json=False,
                                            default_columns=RegulatoryFamilyTransformer.columns)
        apply_processors(list, df=regulatory_family, col='pubmed')
        regulatory_family = regulatory_family.explode('pubmed')

        publication = read_from_stack(tl=self, file='publication', json=False,
                                      default_columns=PublicationTransformer.columns)
        publication = publication.dropna(subset=['pmid'])
        publication = publication.drop_duplicates(subset=['pmid'])

        merged = pd.merge(regulatory_family, publication, left_on='pubmed', right_on='pmid',
                          suffixes=('_regulatory_family', '_publication'))
        merged = merged.dropna(subset=['protrend_id_regulatory_family'])
        merged = merged.dropna(subset=['protrend_id_publication'])
        merged = merged.drop_duplicates(subset=['protrend_id_regulatory_family', 'protrend_id_publication'])

        from_identifiers = merged['protrend_id_regulatory_family'].tolist()
        to_identifiers = merged['protrend_id_publication'].tolist()
        size = len(to_identifiers)

        df = self.make_connection(size=size,
                                  from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)
