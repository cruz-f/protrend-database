import sys

from protrend.runners import Director
from protrend.transform.regprecise import *
from protrend.utils import NeoDatabase, ROOT_PATH

src_path = ROOT_PATH.parent
sys.path.insert(0, str(src_path))
from protrend.runners import run_spider


def transform_runner() -> Director:
    transformers = [
        EffectorTransformer(),
        GeneTransformer(),
        OperonTransformer(),
        OrganismTransformer(),
        PathwayTransformer(),
        PublicationTransformer(),
        RegulatorTransformer(),
        RegulatoryFamilyTransformer(),
        SourceTransformer(),
        TFBSTransformer(),
    ]
    director = Director(transformers=transformers)
    return director


if __name__ == "__main__":
    # run_spider(spider='regprecise', staging_area=STAGING_AREA_PATH, version='0.0.0')

    neo_db = NeoDatabase(user_name='neo4j', password='protrend', ip='localhost', port='7687')
    neo_db.connect()

    transform_director = transform_runner()
    transform_director.transform()
