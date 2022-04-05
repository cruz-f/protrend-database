import pandas as pd
from tqdm import tqdm

from protrend.model import RegulatoryInteraction
from protrend.transform.standardizer.base import StandardizerTransformer


class RegulatoryInteractionTransformer(StandardizerTransformer,
                                       source='standardizer',
                                       version='0.0.0',
                                       node=RegulatoryInteraction,
                                       order=60,
                                       register=True):

    def transform(self):
        nodes = self.fetch_nodes()
        orphans = set()
        for node in tqdm(nodes, desc=f'{self.node.node_name()} - node_standardization'):
            if not node.data_source:
                orphans.add(node.protrend_id)
                node.delete()
                continue

            if node.data_organism and node.data_regulator and node.data_gene:
                continue

            orphans.add(node.protrend_id)
            node.delete()

        df = {'protrend_id': orphans}
        return pd.DataFrame(df)
