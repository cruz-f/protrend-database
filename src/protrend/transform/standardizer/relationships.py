import pandas as pd
from neo4j.exceptions import Neo4jError, DriverError
from tqdm import tqdm

from protrend.model import Organism, Gene, TFBS, RegulatoryInteraction, Regulator
from protrend.transform.standardizer.base import StandardizerTransformer
from protrend.utils import NeoDatabase, Settings


neo_db = NeoDatabase(user_name=Settings.db_user_name, password=Settings.db_password,
                     ip=Settings.db_ip, port=Settings.db_port)


class RelationshipTransformer(StandardizerTransformer,
                              source='standardizer',
                              version='0.0.0',
                              node=Organism,
                              order=100,
                              register=True):

    def transform(self):
        labels = {'Regulator': Regulator, 'Gene': Gene, 'TFBS': TFBS, 'RegulatoryInteraction': RegulatoryInteraction}
        orphans = set()

        for label, node_cls in tqdm(labels.items(), total=len(labels), desc='relationships_standardization'):
            query = f'match(n:{label})-[]-(o:Organism) return n.protrend_id as identifier, count(distinct(o.protrend_id)) as cnt'

            try:
                wrong_connections = neo_db.db.cypher_query(query)

            except (Neo4jError, DriverError):
                wrong_connections = [[], []]

            df = pd.DataFrame(wrong_connections[0], columns=wrong_connections[1])

            if df.empty:
                continue

            df = df[df['cnt'] > 1]

            if df.empty:
                continue

            for identifier in df['identifier']:
                node = node_cls.nodes.get(protrend_id=identifier)

                orphans.add(identifier)
                node.delete()

        df = {'protrend_id': list(orphans)}
        return pd.DataFrame(df)
