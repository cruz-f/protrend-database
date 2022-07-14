import pandas as pd

from protrend.annotation import (PathwayDTO, annotate_pathways)
from protrend.log import ProtrendLogger
from ._utils import get_values


class PathwayMixIn:

    @staticmethod
    def annotate_pathways(df: pd.DataFrame) -> pd.DataFrame:
        input_values = get_values(df, 'input_value')

        pathways = [PathwayDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(pathways)} pathways')

        names = get_values(df, 'name')

        ProtrendLogger.log.info('Annotating with the following params: name')

        annotate_pathways(dtos=pathways, names=names)

        pathways_dict = [dto.to_dict() for dto in pathways]

        pathways_df = pd.DataFrame(pathways_dict)
        return pathways_df
