import pandas as pd
from tqdm import tqdm

from protrend.model import Pathway
from protrend.transform.standardizer.base import StandardizerTransformer


class PathwayTransformer(StandardizerTransformer,
                         source='standardizer',
                         version='0.0.0',
                         node=Pathway,
                         order=50,
                         register=True):

    def transform(self):
        nodes = self.fetch_nodes()
        orphans = set()
        for node in tqdm(nodes, desc=f'{self.node.node_name()} - node_standardization'):
            if not node.data_source:
                orphans.add(node.protrend_id)
                node.delete()
                continue

            if node.regulator or node.gene:
                continue

            orphans.add(node.protrend_id)
            node.delete()

        df = {'protrend_id': list(orphans)}
        return pd.DataFrame(df)
