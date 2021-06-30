from typing import Dict

from neomodel import (StringProperty,
                      IntegerProperty,
                      RelationshipTo,
                      ArrayProperty)

from protrend.models.node import Node
from protrend.models.version import Version


# noinspection PyAbstractClass
class CollecTFVersion(Version):

    database = RelationshipTo('Database', 'VERSIONING')


# noinspection PyAbstractClass
class Database(Node):
    property_as_id = 'name'

    name = StringProperty(default='collectf')
    url = StringProperty(default='http://www.collectf.org/browse/browse/')
    doi = StringProperty(default='10.1093/nar/gkt1123')
    authors = ArrayProperty(StringProperty(), default=['Sefa Kılıç', 'Ivan Erill'])
    description = StringProperty(default='CollecTF: a database of experimentally validated transcription '
                                         'factor-binding sites in Bacteria')

    version = RelationshipTo(CollecTFVersion, 'VERSIONING')


# noinspection PyAbstractClass
class Taxonomy(Node):
    property_as_id = 'taxonomy_id'

    taxonomy_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    organism = RelationshipTo('Organism', 'HAS')


# noinspection PyAbstractClass
class Organism(Node):
    property_as_id = 'name'

    name = StringProperty(required=True)
    genome_accession = StringProperty()

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    regulon = RelationshipTo('Regulon', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class Regulon(Node):
    property_as_id = 'uniprot_accession'

    uniprot_accession = StringProperty(required=True)
    name = StringProperty()
    url = StringProperty()

    organism = RelationshipTo(Organism, 'IS_FROM')
    transcription_factor = RelationshipTo('TranscriptionFactor', 'IS_FROM')
    operon = RelationshipTo('Operon', 'HAS')
    gene = RelationshipTo('Gene', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')
    experimental_evidence = RelationshipTo('ExperimentalEvidence', 'IS_FROM')


# noinspection PyAbstractClass
class Operon(Node):
    property_as_id = 'operon_id'

    operon_id = StringProperty(required=True)

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    gene = RelationshipTo('Gene', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class Gene(Node):
    property_as_id = 'locus_tag'

    locus_tag = StringProperty(required=True)

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class TFBS(Node):
    property_as_id = 'tfbs_id'

    tfbs_id = StringProperty(required=True)
    site_start = IntegerProperty(required=True)
    site_end = IntegerProperty(required=True)
    site_strand = IntegerProperty(required=True)
    sequence = StringProperty(required=True)
    mode = StringProperty(required=True)
    pubmed = ArrayProperty(StringProperty())

    organism = RelationshipTo(Organism, 'IS_FROM')
    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    gene = RelationshipTo(Gene, 'IS_FROM')
    experimental_evidence = RelationshipTo('ExperimentalEvidence', 'IS_FROM')


# noinspection PyAbstractClass
class ExperimentalEvidence(Node):
    property_as_id = 'exp_id'

    exp_id = StringProperty(required=True)


# noinspection PyAbstractClass
class TranscriptionFactor(Node):
    property_as_id = 'name'

    name = StringProperty(required=True)
    family = StringProperty()
    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    regulon = RelationshipTo(Regulon, 'HAS')
