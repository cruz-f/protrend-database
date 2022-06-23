import pandas as pd
import requests
from neo4j.exceptions import Neo4jError, DriverError

from .base import FunctionalTFBSTransformer
from protrend.model import TFBS
from ...utils import Settings


class TFBSTransformer(FunctionalTFBSTransformer,
                      source='functional_tfbs',
                      version='0.0.0',
                      node=TFBS,
                      order=90,
                      register=True):

    def fetch_nodes(self):
        try:
            regulators = Regulator.nodes.all()

            identifiers = []
            loci = []
            names = []
            binding_sites = []
            for regulator in regulators:
                sites = regulator.tfbs.all()
                if not sites:
                    continue
                identifiers.append(regulator.protrend_id)
                loci.append(regulator.locus_tag)
                names.append(regulator.name)
                binding_sites.append([site.sequence for site in sites])

            df = pd.DataFrame(data={'protrend_id': identifiers,
                                    'locus_tag': loci,
                                    'name': names,
                                    'binding_sites': binding_sites})
            return df
        except (Neo4jError, DriverError) as e:
            print(e)
            return pd.DataFrame()

    def align_tfbs(self, sequences: list, headers: list) -> pd.DataFrame:

        data = {
            'sequences': sequences,
            'headers': headers
        }

        payload = {'sequences': data, 'k': 0}

        response = requests.request("POST", Settings.lasagna_url, json=payload)

        return pd.Dataframe(response.json())

    def transform(self) -> pd.DataFrame:
        tfbs = self.fetch_nodes()

        sequences = tfbs["tfbs"].to_list()
        headers = tfbs["regulator"].to_list()

        final_tfbs = []

        for seqs, heads in zip(sequences, headers):
            aligned_tfbs = self.align_tfbs(seqs, heads)
            final_tfbs.append(aligned_tfbs)

        self.stack_transformed_nodes(final_tfbs)
        return final_tfbs