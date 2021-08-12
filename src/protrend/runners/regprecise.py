import sys

from protrend.runners.director import Director
from protrend.transform.regprecise import *
from protrend.utils.db_connection import DBSettings

sys.path.insert(0, r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src')


# from protrend.runners.spider_runner import run_spider

def transform_runner() -> Director:
    transformers = [
        # EffectorTransformer(),
        GeneTransformer(),
        OperonTransformer(),
        # OrganismTransformer(),
        # PathwayTransformer(),
        # PublicationTransformer(),
        RegulatorTransformer(),
        # RegulatoryFamilyTransformer(),
        # SourceTransformer(),
        TFBSTransformer(),
    ]
    director = Director(transformers=transformers)
    return director


if __name__ == "__main__":
    # run_spider(spider='regprecise',
    #            staging_area=r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src\protrend\extract\staging_area',
    #            version='0.0.0')

    db_settings = DBSettings(user_name='neo4j',
                             password='protrend',
                             ip='localhost',
                             port='7687',
                             db_name='neo4j',
                             dbms='')
    db_settings.connect()

    transform_director = transform_runner()
    transform_director.transform()
