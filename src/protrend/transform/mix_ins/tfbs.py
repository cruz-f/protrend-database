import pandas as pd

from protrend.utils import apply_processors
from protrend.utils.processors import to_list_nan, protrend_hash
from protrend.transform.transformations import drop_duplicates, drop_empty_string


class TFBSMixIn:

    @staticmethod
    def site_hash(df: pd.DataFrame) -> pd.DataFrame:
        """
        It creates a site hash using organism + sequence + strand + start + stop + length
        """

        # filter by organism + sequence + strand + start + stop + length
        df2 = apply_processors(df,
                               organism=to_list_nan,
                               sequence=to_list_nan,
                               strand=to_list_nan,
                               start=to_list_nan,
                               stop=to_list_nan,
                               length=to_list_nan)

        ri_series_hash = df2['organism'].copy()
        ri_series_hash += df2['sequence'].copy()
        ri_series_hash += df2['strand'].copy()
        ri_series_hash += df2['start'].copy()
        ri_series_hash += df2['stop'].copy()
        ri_series_hash += df2['length'].copy()

        df = df.assign(site_hash=ri_series_hash)
        df = apply_processors(df, site_hash=protrend_hash)
        df = df.dropna(subset=['site_hash'])
        df = drop_empty_string(df, 'site_hash')
        df = drop_duplicates(df=df, subset=['site_hash'])

        return df
