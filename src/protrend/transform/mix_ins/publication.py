import pandas as pd

from protrend.annotation import (PublicationDTO, annotate_publications)
from protrend.log import ProtrendLogger
from ._utils import get_values


class PublicationMixIn:

    @staticmethod
    def annotate_publications(df: pd.DataFrame) -> pd.DataFrame:

        input_values = get_values(df, 'input_value')

        publications = [PublicationDTO(input_value=input_value) for input_value in input_values]

        ProtrendLogger.log.info(f'Annotating {len(publications)} publications')

        identifiers = get_values(df, 'pmid')

        ProtrendLogger.log.info('Annotating with the following params: pmid')

        annotate_publications(dtos=publications, identifiers=identifiers)

        publications_dict = [dto.to_dict() for dto in publications]

        publications_df = pd.DataFrame(publications_dict)
        return publications_df
