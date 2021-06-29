from typing import Dict

from neomodel import (StringProperty,
                      IntegerProperty,
                      RelationshipTo,
                      ArrayProperty,
                      FloatProperty)

from protrend.models.node import Node
from protrend.models.version import Version


# noinspection PyAbstractClass
class RegPreciseVersion(Version):

    database = RelationshipTo('Database', 'VERSIONING')


# noinspection PyAbstractClass
class Database(Node):
    property_as_id = 'name'

    name = StringProperty(default='regprecise')
    url = StringProperty(default='https://regprecise.lbl.gov/collections.jsp')
    doi = StringProperty(default='10.1186/1471-2164-14-745')
    authors = ArrayProperty(StringProperty(), default=['Pavel S Novichkov', 'Dmitry A Rodionov'])
    description = StringProperty(default='Collection of Manually Curated Inferences of Regulons in Prokaryotic Genomes')

    version = RelationshipTo(RegPreciseVersion, 'VERSIONING')


# noinspection PyAbstractClass
class Taxonomy(Node):
    property_as_id = 'collection_id'

    collection_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    genome = RelationshipTo('Genome', 'HAS')


# noinspection PyAbstractClass
class Genome(Node):
    property_as_id = 'genome_id'

    genome_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    regulon = RelationshipTo('Regulon', 'HAS')


# noinspection PyAbstractClass
class TranscriptionFactor(Node):
    property_as_id = 'collection_id'

    collection_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Regulog(Node):
    property_as_id = 'regulog_id'

    regulog_id = IntegerProperty(required=True)
    name = StringProperty()
    url = StringProperty()

    regulator_type = StringProperty()
    regulator_family = StringProperty()
    rfam = StringProperty()
    regulation_mode = StringProperty()
    biological_process = StringProperty()
    regulation_effector = StringProperty()
    phylum = StringProperty()

    regulon = RelationshipTo('Regulon', 'HAS')

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    transcription_factor = RelationshipTo(TranscriptionFactor, 'IS_FROM')
    tf_family = RelationshipTo('TranscriptionFactorFamily', 'IS_FROM')
    rna_family = RelationshipTo('RNAFamily', 'IS_FROM')
    effector = RelationshipTo('Effector', 'IS_FROM')
    pathway = RelationshipTo('Pathway', 'IS_FROM')


# noinspection PyAbstractClass
class TranscriptionFactorFamily(Node):
    property_as_id = 'tffamily_id'

    tffamily_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class RNAFamily(Node):
    property_as_id = 'riboswitch_id'

    riboswitch_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())
    rfam = StringProperty()

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Effector(Node):
    property_as_id = 'effector_id'

    effector_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Pathway(Node):
    property_as_id = 'pathway_id'

    pathway_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Regulon(Node):
    property_as_id = 'regulon_id'

    regulon_id = IntegerProperty(required=True)
    name = StringProperty()
    url = StringProperty()

    regulator_type = StringProperty()
    regulator_locus_tag = StringProperty()
    rfam = StringProperty()
    regulator_family = StringProperty()
    regulation_mode = StringProperty()
    biological_process = StringProperty()
    regulation_effector = StringProperty()
    regulation_regulog = StringProperty()

    genome = RelationshipTo(Genome, 'IS_FROM')
    operon = RelationshipTo('Operon', 'HAS')
    gene = RelationshipTo('Gene', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')

    regulog = RelationshipTo(Regulog, 'IS_FROM')
    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    transcription_factor = RelationshipTo(TranscriptionFactor, 'IS_FROM')
    tf_family = RelationshipTo(TranscriptionFactorFamily, 'IS_FROM')
    rna_family = RelationshipTo(RNAFamily, 'IS_FROM')
    effector = RelationshipTo(Effector, 'IS_FROM')
    pathway = RelationshipTo(Pathway, 'IS_FROM')


# noinspection PyAbstractClass
class Operon(Node):
    property_as_id = 'operon_id'

    operon_id = StringProperty(required=True)
    name = StringProperty()
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    gene = RelationshipTo('Gene', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class Gene(Node):
    property_as_id = 'locus_tag'

    locus_tag = StringProperty(required=True)
    name = StringProperty()
    function = StringProperty()
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class TFBS(Node):
    property_as_id = 'tfbs_id'

    tfbs_id = StringProperty(required=True)
    position = IntegerProperty(required=True)
    score = FloatProperty(required=True)
    sequence = StringProperty(required=True)
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    gene = RelationshipTo(Gene, 'IS_FROM')
