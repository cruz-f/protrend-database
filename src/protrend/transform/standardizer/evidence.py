import pandas as pd
from tqdm import tqdm

from protrend.model import Evidence
from protrend.transform.standardizer.base import StandardizerTransformer


class EvidenceTransformer(StandardizerTransformer,
                          source='standardizer',
                          version='0.0.0',
                          node=Evidence,
                          order=40,
                          register=True):

    def transform(self):
        nodes = self.fetch_nodes()
        orphans = set()
        for node in tqdm(nodes, desc=f'{self.node.node_name()} - node_standardization'):
            if node.regulator or \
                    node.operon or \
                    node.gene or \
                    node.tfbs or \
                    node.regulatory_interaction:
                continue

            orphans.add(node.protrend_id)
            node.delete()

        df = {'protrend_id': list(orphans)}
        return pd.DataFrame(df)
