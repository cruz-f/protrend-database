import pandas as pd

from protrend.annotation import (EffectorDTO, annotate_effectors)
from protrend.log import ProtrendLogger
from ._utils import get_values


class EffectorMixIn:

    @staticmethod
    def annotate_effectors(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        effectors = [EffectorDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(effectors)} effectors')

        names = get_values(df, 'name')

        ProtrendLogger.log.info('Annotating with the following params: name')

        annotate_effectors(dtos=effectors, names=names)

        effectors_dict = [dto.to_dict() for dto in effectors]

        effectors_df = pd.DataFrame(effectors_dict)
        return effectors_df
