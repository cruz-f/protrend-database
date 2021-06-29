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

    taxonomy = RelationshipTo('Taxonomy', 'VERSIONING')
    genome = RelationshipTo('Genome', 'VERSIONING')

    transcriptionfactor = RelationshipTo('TranscriptionFactor', 'VERSIONING')
    regulog = RelationshipTo('Regulog', 'VERSIONING')

    transcriptionfactorfamily = RelationshipTo('TranscriptionFactorFamily', 'VERSIONING')
    rnafamily = RelationshipTo('RNAFamily', 'VERSIONING')
    effector = RelationshipTo('Effector', 'VERSIONING')
    pathway = RelationshipTo('Pathway', 'VERSIONING')

    regulon = RelationshipTo('Regulon', 'VERSIONING')
    operon = RelationshipTo('Operon', 'VERSIONING')
    gene = RelationshipTo('Gene', 'VERSIONING')
    tfbs = RelationshipTo('TFBS', 'VERSIONING')

    @property
    def versioned_nodes(self) -> Dict[str, RelationshipTo]:

        return {'Database': self.database,
                'Taxonomy': self.taxonomy,
                'Genome': self.genome,
                'TranscriptionFactor': self.transcriptionfactor,
                'Regulog': self.regulog,
                'TranscriptionFactorFamily': self.transcriptionfactorfamily,
                'RNAFamily': self.rnafamily,
                'Effector': self.effector,
                'Pathway': self.pathway,
                'Regulon': self.regulon,
                'Operon': self.operon,
                'Gene': self.gene,
                'TFBS': self.tfbs}


class RegPreciseNode(Node):

    __abstract_node__ = True

    version = RelationshipTo(RegPreciseVersion, 'VERSIONING')


# noinspection PyAbstractClass
class Database(RegPreciseNode):
    property_as_id = 'name'

    name = StringProperty(default='regprecise')
    url = StringProperty(default='https://regprecise.lbl.gov/collections.jsp')
    doi = StringProperty(default='10.1186/1471-2164-14-745')
    authors = ArrayProperty(StringProperty(), default=['Pavel S Novichkov', 'Dmitry A Rodionov'])
    description = StringProperty(default='Collection of Manually Curated Inferences of Regulons in Prokaryotic Genomes')


# noinspection PyAbstractClass
class Taxonomy(RegPreciseNode):
    property_as_id = 'collection_id'

    collection_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    genome = RelationshipTo('Genome', 'HAS')


# noinspection PyAbstractClass
class Genome(RegPreciseNode):
    property_as_id = 'genome_id'

    genome_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    regulon = RelationshipTo('Regulon', 'HAS')


# noinspection PyAbstractClass
class TranscriptionFactor(RegPreciseNode):
    property_as_id = 'collection_id'

    collection_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Regulog(RegPreciseNode):
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
class TranscriptionFactorFamily(RegPreciseNode):
    property_as_id = 'tffamily_id'

    tffamily_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class RNAFamily(RegPreciseNode):
    property_as_id = 'riboswitch_id'

    riboswitch_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())
    rfam = StringProperty()

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Effector(RegPreciseNode):
    property_as_id = 'effector_id'

    effector_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Pathway(RegPreciseNode):
    property_as_id = 'pathway_id'

    pathway_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulog = RelationshipTo('Regulog', 'HAS')


# noinspection PyAbstractClass
class Regulon(RegPreciseNode):
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
class Operon(RegPreciseNode):
    property_as_id = 'operon_id'

    operon_id = StringProperty(required=True)
    name = StringProperty()
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    gene = RelationshipTo('Gene', 'HAS')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class Gene(RegPreciseNode):
    property_as_id = 'locus_tag'

    locus_tag = StringProperty(required=True)
    name = StringProperty()
    function = StringProperty()
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    tfbs = RelationshipTo('TFBS', 'HAS')


# noinspection PyAbstractClass
class TFBS(RegPreciseNode):
    property_as_id = 'tfbs_id'

    tfbs_id = StringProperty(required=True)
    position = IntegerProperty(required=True)
    score = FloatProperty(required=True)
    sequence = StringProperty(required=True)
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    operon = RelationshipTo(Operon, 'IS_FROM')
    gene = RelationshipTo(Gene, 'IS_FROM')
