import pandas as pd
import requests
from neo4j.exceptions import Neo4jError, DriverError

from .base import FunctionalTFBSTransformer
from protrend.model import TFBS
from protrend.io.utils import read_promoters


class TFBSTransformer(FunctionalTFBSTransformer,
                      source='functional_tfbs',
                      version='0.0.0',
                      node=TFBS,
                      order=90,
                      register=True):

    def fetch_nodes(self):
        try:
            return self.node.node_to_df()
        except (Neo4jError, DriverError) as e:
            print(e)
            return pd.DataFrame()

    def align_tfbs(self, tfbs: pd.DataFrame, promoters: pd.DataFrame) -> pd.DataFrame:
        url = "http://127.0.0.1:5000/submit"

        sequences = {
            'sequences': [
                'GACCAATGATTATTAGCCAAT',
                'GGCCAGCCTTGCCTTGACCAATAGCCTTGACAAGGCAAACTT',
                'CCCGAGCCGCTGATTGGCTTTCAGG',
                'ACCAATGACTTTTAAGTACC',
                'CCAAT'

            ],
            'headers': [
                'site_0',
                'site_1',
                'site_2',
                'site_3',
                'site_4'
            ]}

        payload = {'sequences': sequences, 'k': 0}

        response = requests.request("POST", url, json=payload)

        #print(response.json())
        return response.json()

    def calculate_descriptors(self, aligned_tfbs: pd.DataFrame) -> pd.DataFrame:
        # TODO: method to calculate descriptors - pwm, pssm, gc content
        pass

    def transform(self) -> pd.DataFrame:
        tfbs = self.fetch_nodes()
        promoters = read_promoters(source=self.source, version=self.version, columns=[])

        aligned_tfbs = self.align_tfbs(tfbs, promoters)

        final_tfbs = self.calculate_descriptors(aligned_tfbs)

        self.stack_transformed_nodes(final_tfbs)
        return final_tfbs
