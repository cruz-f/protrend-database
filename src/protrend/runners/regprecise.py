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

    # TODO: Regulator:
    #  - wrong name in some regulators
    #  - parse uniprot query misses some valid uniprot entries - CORRECTED
    #  - strand is missing
    #  - synonyms are not unique - CORRECTED
    #  - some regulators are not being integrated and being dropped

    # TODO: Gene:
    #  - wrong name in some genes
    #  - parse uniprot query misses some valid uniprot entries - CORRECTED
    #  - strand is missing
    #  - synonyms are not unique - CORRECTED
    #  - missing most UniProt accessions in gene annotation
    #  - regulon column is wrongly parsed
    #  - missing all information regarding organism and ncbi taxonomy
    #  - some genes are not being integrated and being dropped

    # TODO: Operon:
    #  - Missing all hits with tfbs (might not be a problem due to sampling)
    #  - All genes were found (might not be a problem due to sampling)
    #  - regulon column is wrongly parsed as set of list of str
    #  - first position left is wrongly inferred
    #  - genes_id column should be dropped
    #  - strand is missing

    # TODO: RegulatoryFamily:
    #  - Description is missing in some rows

    # TODO: Connectors
