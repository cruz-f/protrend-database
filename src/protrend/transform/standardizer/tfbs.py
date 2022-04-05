import pandas as pd
from tqdm import tqdm

from protrend.model import TFBS
from protrend.transform.standardizer.base import StandardizerTransformer


class TFBSTransformer(StandardizerTransformer,
                      source='standardizer',
                      version='0.0.0',
                      node=TFBS,
                      order=70,
                      register=True):

    def transform(self):
        nodes = self.fetch_nodes()
        orphans = set()
        for node in tqdm(nodes, desc=f'{self.node.node_name()} - node_standardization'):
            if not node.data_source:
                orphans.add(node.protrend_id)
                node.delete()
                continue

            if node.data_organism and node.regulator and node.gene and node.regulatory_interaction:
                continue

            orphans.add(node.protrend_id)
            node.delete()

        df = {'protrend_id': orphans}
        return pd.DataFrame(df)
