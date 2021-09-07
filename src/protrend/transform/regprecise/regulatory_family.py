import pandas as pd

from protrend.io.json import read_json_lines, read_json_frame
from protrend.io.utils import read_from_stack
from protrend.transform.connector import Connector
from protrend.transform.processors import (remove_white_space, remove_regprecise_more, remove_multiple_white_space,
                                           rstrip, lstrip, remove_pubmed, apply_processors, to_set, to_list_nan,
                                           to_int_str, to_nan, to_list)
from protrend.transform.regprecise import PublicationTransformer
from protrend.transform.regprecise.regulator import RegulatorTransformer
from protrend.transform.regprecise.settings import (RegulatoryFamilySettings, RegulatoryFamilyToSource,
                                                    RegulatoryFamilyToPublication, RegulatoryFamilyToRegulator)
from protrend.transform.regprecise.source import SourceTransformer
from protrend.transform.transformer import Transformer


class RegulatoryFamilyTransformer(Transformer):
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

        df = apply_processors(df, name=remove_white_space,
                              description=[remove_regprecise_more, remove_pubmed, remove_multiple_white_space, rstrip,
                                           lstrip])

        return df

    def _transform_tf(self, tf: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=tf, subset=['name'], perfect_match=True, preserve_nan=False)

        df = apply_processors(df, name=remove_white_space,
                              description=[remove_regprecise_more, remove_pubmed, remove_multiple_white_space, rstrip,
                                           lstrip])

        return df

    def _transform_rna(self, rna: pd.DataFrame) -> pd.DataFrame:
        df = self.drop_duplicates(df=rna, subset=['rfam'], perfect_match=True, preserve_nan=False)

        df = apply_processors(df, name=remove_white_space, rfam=remove_white_space,
                              description=[remove_regprecise_more, remove_pubmed, remove_multiple_white_space, rstrip,
                                           lstrip],
                              pubmed=to_set,
                              regulog=to_set)

        df = df.rename(columns={'url': 'url_rna'})

        # set mechanism
        df['mechanism'] = 'small RNA (sRNA)'
        return df

    def _transform_tfs(self, tf_family: pd.DataFrame, tf: pd.DataFrame) -> pd.DataFrame:
        df = pd.merge(tf_family, tf, how='outer', on='name', suffixes=('_tf_family', '_tf'))

        # set mechanism
        df['mechanism'] = 'transcription factor'

        # concat description
        df = apply_processors(df, description_tf_family=to_nan, description_tf=to_nan)
        df['description_tf_family'] = df['description_tf_family'].fillna(value='')
        df['description_tf'] = df['description_tf'].fillna(value='')
        df = self.concat_columns(df=df, column='description', left='description_tf_family', right='description_tf')

        # concat pubmed
        df = apply_processors(df, pubmed_tf_family=to_list_nan, pubmed_tf=to_list_nan)
        df = self.concat_columns(df=df, column='pubmed', left='pubmed_tf_family', right='pubmed_tf')
        df = apply_processors(df, pubmed=to_set)

        # concat regulog
        df = apply_processors(df, regulog_tf_family=to_list_nan, regulog_tf=to_list_nan)
        df = self.concat_columns(df=df, column='regulog', left='regulog_tf_family', right='regulog_tf')
        df = apply_processors(df, regulog=to_set)

        return df

    def transform(self):
        # -------------------- TFs -------------------
        tf_family = read_from_stack(stack=self._transform_stack, file='tf_family', 
                                    default_columns=self.tf_family_columns, reader=read_json_lines)
        tf_family = self._transform_tf_family(tf_family)

        tf = read_from_stack(stack=self._transform_stack, file='tf', 
                             default_columns=self.tf_columns, reader=read_json_lines)
        tf = self._transform_tf(tf)
        tfs = self._transform_tfs(tf_family=tf_family, tf=tf)

        # -------------------- RNAs -------------------
        rna = read_from_stack(stack=self._transform_stack, file='rna', 
                              default_columns=self.rna_columns, reader=read_json_lines)
        rna = self._transform_rna(rna)

        df = pd.concat([tfs, rna])

        df = apply_processors(df, tffamily_id=to_int_str, riboswitch_id=to_int_str, collection_id=to_int_str)

        self._stack_transformed_nodes(df)
        return df


class RegulatoryFamilyToSourceConnector(Connector):
    default_settings = RegulatoryFamilyToSource

    def connect(self):
        regulatory_family = read_from_stack(stack=self._connect_stack, file='regulatory_family',
                                            default_columns=RegulatoryFamilyTransformer.columns,
                                            reader=read_json_frame)
        source = read_from_stack(stack=self._connect_stack, file='source',
                                 default_columns=SourceTransformer.columns, reader=read_json_frame)

        protrend_id = source['protrend_id'].iloc[0]

        tf_chain = [('tffamily_id', 'url_tf_family'),
                    ('collection_id', 'url_tf'),
                    ('riboswitch_id', 'url_rna')]

        dfs = []
        for key, url in tf_chain:
            key_df = regulatory_family.dropna(subset=[key])

            from_identifiers = key_df['protrend_id'].tolist()
            size = len(from_identifiers)

            to_identifiers = [protrend_id] * size

            kwargs = dict(url=key_df[url].tolist(),
                          external_identifier=key_df[key].tolist(),
                          key=[key] * size)

            df = self.make_connection(from_identifiers=from_identifiers,
                                      to_identifiers=to_identifiers,
                                      kwargs=kwargs)

            dfs.append(df)

        df = pd.concat(dfs, axis=0)

        self.stack_csv(df)


class RegulatoryFamilyToPublicationConnector(Connector):
    default_settings = RegulatoryFamilyToPublication

    def connect(self):
        regulatory_family = read_from_stack(stack=self._connect_stack, file='regulatory_family',
                                            default_columns=RegulatoryFamilyTransformer.columns, reader=read_json_frame)
        regulatory_family = apply_processors(regulatory_family, pubmed=to_list)
        regulatory_family = regulatory_family.explode('pubmed')

        publication = read_from_stack(stack=self._connect_stack, file='publication',
                                      default_columns=PublicationTransformer.columns, reader=read_json_frame)
        publication = publication.dropna(subset=['pmid'])
        publication = publication.drop_duplicates(subset=['pmid'])

        merged = pd.merge(regulatory_family, publication, left_on='pubmed', right_on='pmid',
                          suffixes=('_regulatory_family', '_publication'))
        merged = merged.dropna(subset=['protrend_id_regulatory_family'])
        merged = merged.dropna(subset=['protrend_id_publication'])
        merged = merged.drop_duplicates(subset=['protrend_id_regulatory_family', 'protrend_id_publication'])

        from_identifiers = merged['protrend_id_regulatory_family'].tolist()
        to_identifiers = merged['protrend_id_publication'].tolist()

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)


class RegulatoryFamilyToRegulatorConnector(Connector):
    default_settings = RegulatoryFamilyToRegulator

    def connect(self):
        regulatory_family = read_from_stack(stack=self._connect_stack, file='regulatory_family',
                                            default_columns=RegulatoryFamilyTransformer.columns, reader=read_json_frame)
        regulator = read_from_stack(stack=self._connect_stack, file='regulator',
                                    default_columns=RegulatorTransformer.columns, reader=read_json_frame)

        tfs_chain = [('tffamily_id', 'tf_family'),
                     ('collection_id', 'transcription_factor'),
                     ('riboswitch_id', 'rna_family')]

        from_identifiers = []
        to_identifiers = []

        for family_key, regulator_key in tfs_chain:
            family_df = regulatory_family.dropna(subset=[family_key])

            regulator_df = regulator.explode(regulator_key)
            regulator_df = regulator_df.dropna(subset=[regulator_key])

            merged = pd.merge(family_df, regulator_df, left_on=family_key, right_on=regulator_key,
                              suffixes=('_regulatory_family', '_regulator'))

            merged = merged.dropna(subset=['protrend_id_regulatory_family'])
            merged = merged.dropna(subset=['protrend_id_regulator'])
            merged = merged.drop_duplicates(subset=['protrend_id_regulatory_family', 'protrend_id_regulator'])

            from_identifiers.extend(merged['protrend_id_regulatory_family'].tolist())
            to_identifiers.extend(merged['protrend_id_regulator'].tolist())

        df = self.make_connection(from_identifiers=from_identifiers,
                                  to_identifiers=to_identifiers)

        self.stack_csv(df)
