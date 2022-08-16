import pandas as pd
from tqdm import tqdm

from protrend.model import Gene
from protrend.transform.standardizer.base import StandardizerTransformer


class GeneTransformer(StandardizerTransformer,
                      source='standardizer',
                      version='0.0.0',
                      node=Gene,
                      order=80,
                      register=True):

    def transform(self):
        nodes = self.fetch_nodes()
        orphans = set()
        for node in tqdm(nodes, desc=f'{self.node.node_name()} - node_standardization'):
            if not node.data_source:
                orphans.add(node.protrend_id)
                node.delete()
                continue

            if node.operon:
                continue

            if node.organism and node.regulator and node.regulatory_interaction:
                continue

            orphans.add(node.protrend_id)
            node.delete()

        df = {'protrend_id': list(orphans)}
        return pd.DataFrame(df)
