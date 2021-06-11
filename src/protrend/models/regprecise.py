from neomodel import (StructuredNode,
                      StringProperty,
                      IntegerProperty,
                      RelationshipTo,
                      ArrayProperty,
                      FloatProperty)

from protrend.models.node import Node
from protrend.models.version import VersionNode


# noinspection PyAbstractClass
class Version(StructuredNode, VersionNode):
    pass


# noinspection PyAbstractClass
class Collection(StructuredNode, Node):
    name = StringProperty(default='regprecise collection')
    url = StringProperty(default='https://regprecise.lbl.gov/collections.jsp')
    doi = StringProperty(default='10.1186/1471-2164-14-745')
    authors = ArrayProperty(StringProperty(), default=['Pavel S Novichkov', 'Dmitry A Rodionov'])
    description = StringProperty(default='Collection of Manually Curated Inferences of Regulons in Prokaryotic Genomes')

    taxonomy = RelationshipTo('Taxonomy', 'HAS')
    transcription_factor = RelationshipTo('TranscriptionFactor', 'HAS')
    transcription_factor_family = RelationshipTo('TranscriptionFactorFamily', 'HAS')
    rna_family = RelationshipTo('RNAFamily', 'HAS')
    effector = RelationshipTo('Effector', 'HAS')
    pathway = RelationshipTo('Pathway', 'HAS')


# noinspection PyAbstractClass
class Taxonomy(StructuredNode, Node):
    collection_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    collection = RelationshipTo(Collection, 'IS_FROM')
    genome = RelationshipTo('Genome', 'HAS')


# noinspection PyAbstractClass
class Genome(StructuredNode, Node):
    genome_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    regulon = RelationshipTo('Regulon', 'HAS')


# noinspection PyAbstractClass
class TranscriptionFactor(StructuredNode, Node):
    collection_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    collection = RelationshipTo(Collection, 'IS_FROM')
    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Regulog(StructuredNode, Node):
    regulog_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulator_type = StringProperty()
    regulator_family = StringProperty()
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
class TranscriptionFactorFamily(StructuredNode, Node):
    tffamily_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    collection = RelationshipTo(Collection, 'IS_FROM')
    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class RNAFamily(StructuredNode, Node):
    riboswitch_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())
    rfam = StringProperty()

    collection = RelationshipTo(Collection, 'IS_FROM')
    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Effector(StructuredNode, Node):
    effector_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    collection = RelationshipTo(Collection, 'IS_FROM')
    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Pathway(StructuredNode, Node):
    pathway_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    collection = RelationshipTo(Collection, 'IS_FROM')
    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Regulon(StructuredNode, Node):
    regulon_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulator_type = StringProperty()
    regulator_locus_tag = StringProperty(required=True)
    regulator_family = StringProperty()
    regulation_mode = StringProperty()
    biological_process = StringProperty()
    regulation_effector = StringProperty()
    regulation_regulog = StringProperty(required=True)

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
class Operon(StructuredNode, Node):
    operon_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    gene = RelationshipTo('Gene', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class Gene(StructuredNode, Node):
    locus_tag = StringProperty(required=True)
    name = StringProperty()
    function = StringProperty()
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class TFBS(StructuredNode, Node):
    tfbs_id = StringProperty(required=True)
    position = IntegerProperty(required=True)
    score = FloatProperty(required=True)
    sequence = StringProperty(required=True)
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    gene = RelationshipTo(Gene, 'IS_FROM')
