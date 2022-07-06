from abc import abstractmethod
from typing import Union, Callable

import pandas as pd

from protrend.utils import apply_processors
from protrend.utils.processors import to_list_nan, protrend_hash
from protrend.transform.transformations import drop_duplicates
from protrend.transform.transformer import Transformer


class RegulatoryInteractionMixIn:

    @staticmethod
    def interaction_hash(df: pd.DataFrame) -> pd.DataFrame:
        # filter by organism + regulator + gene + tfbs + effector + regulatory effect
        df2 = apply_processors(df,
                               organism=to_list_nan,
                               regulator=to_list_nan,
                               gene=to_list_nan,
                               tfbs=to_list_nan,
                               effector=to_list_nan,
                               regulatory_effect=to_list_nan)

        ri_series_hash = df2['organism'].copy()
        ri_series_hash += df2['regulator'].copy()
        ri_series_hash += df2['gene'].copy()
        ri_series_hash += df2['tfbs'].copy()
        ri_series_hash += df2['effector'].copy()
        ri_series_hash += df2['regulatory_effect'].copy()

        df = df.assign(interaction_hash=ri_series_hash)
        df = apply_processors(df, interaction_hash=protrend_hash)
        df = drop_duplicates(df=df, subset=['interaction_hash'], perfect_match=True)
        df = df.dropna(subset=['interaction_hash'])

        return df

    @abstractmethod
    def transform_network(self, network: pd.DataFrame) -> pd.DataFrame:
        pass

    @abstractmethod
    def transform_organism(self, organism: pd.DataFrame) -> pd.DataFrame:
        pass

    @abstractmethod
    def transform_regulator(self, regulator: pd.DataFrame) -> pd.DataFrame:
        pass

    @abstractmethod
    def transform_gene(self, gene: pd.DataFrame) -> pd.DataFrame:
        pass

    def transform_tfbs(self, tfbs: pd.DataFrame) -> pd.DataFrame:
        pass

    def transform_effector(self, effector: pd.DataFrame) -> pd.DataFrame:
        pass

    def _transform(self: Union[Transformer, 'RegulatoryInteractionMixIn'],
                   network: pd.DataFrame,
                   organism: pd.DataFrame,
                   organism_key: str,
                   regulator: pd.DataFrame,
                   regulator_key: str,
                   gene: pd.DataFrame,
                   gene_key: str,
                   tfbs: pd.DataFrame = None,
                   tfbs_key: str = None,
                   effector: pd.DataFrame = None,
                   effector_key: str = None,
                   regulatory_effect_processor: Callable = None) -> pd.DataFrame:
        network = self.transform_network(network)

        organism = self.transform_organism(organism)
        regulator = self.transform_regulator(regulator)
        gene = self.transform_gene(gene)

        if tfbs is not None:
            tfbs = self.transform_tfbs(tfbs)

        if effector is not None:
            effector = self.transform_effector(effector)

        regulatory_interaction = pd.merge(network, organism, on=organism_key)
        regulatory_interaction = pd.merge(regulatory_interaction, regulator, on=regulator_key)
        regulatory_interaction = pd.merge(regulatory_interaction, gene, on=gene_key)

        if tfbs is not None:
            regulatory_interaction = pd.merge(regulatory_interaction, tfbs, how='left', on=tfbs_key)

        else:
            regulatory_interaction = regulatory_interaction.assign(tfbs=None)

        if effector is not None:
            regulatory_interaction = pd.merge(regulatory_interaction, effector, how='left', on=effector_key)

        else:
            regulatory_interaction = regulatory_interaction.assign(effector=None)

        if regulatory_effect_processor is not None:
            regulatory_interaction = apply_processors(regulatory_interaction,
                                                      regulatory_effect=regulatory_effect_processor)

        regulatory_interaction = self.interaction_hash(regulatory_interaction)
        return regulatory_interaction
