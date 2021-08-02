import sys

from protrend.runners.director import Director
from protrend.transform.regprecise.organism import OrganismTransformer
from protrend.transform.regprecise.settings import OrganismSettings
from protrend.utils.db_connection import DBSettings

sys.path.insert(0, r'C:\Users\BiSBII\OneDrive - Universidade do Minho\PhD\Protrend\main\protrend-database\src')


# from protrend.runners.spider_runner import run_spider


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

    organisms_settings = OrganismSettings(files={'genome': 'Genome-test.json'})
    transformer = OrganismTransformer(organisms_settings)
    director = Director(transformers=[transformer])
    director.transform()
    print()
