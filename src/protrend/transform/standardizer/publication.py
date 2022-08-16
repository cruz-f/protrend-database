import pandas as pd
from tqdm import tqdm

from protrend.model import Publication
from protrend.transform.standardizer.base import StandardizerTransformer


class PublicationTransformer(StandardizerTransformer,
                             source='standardizer',
                             version='0.0.0',
                             node=Publication,
                             order=40,
                             register=True):

    def transform(self):
        nodes = self.fetch_nodes()
        orphans = set()
        for node in tqdm(nodes, desc=f'{self.node.node_name()} - node_standardization'):
            if node.regulatory_family or \
                    node.regulator or \
                    node.operon or \
                    node.gene or \
                    node.tfbs or \
                    node.regulatory_interaction:
                continue

            orphans.add(node.protrend_id)
            node.delete()

        df = {'protrend_id': list(orphans)}
        return pd.DataFrame(df)
