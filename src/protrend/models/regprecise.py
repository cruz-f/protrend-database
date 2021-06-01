from typing import Union

from neomodel import (StructuredNode,
                      StringProperty,
                      IntegerProperty,
                      RelationshipTo,
                      ArrayProperty,
                      FloatProperty,
                      OneOrMore,
                      UniqueIdProperty,
                      DateTimeProperty, RelationshipFrom)


# noinspection PyAbstractClass
class Version(StructuredNode):
    uid = UniqueIdProperty()
    name = StringProperty(required=True, unique_index=True)
    created = DateTimeProperty(default_now=True)

    details = RelationshipFrom('Details', 'VERSIONING')
    collection = RelationshipFrom('Collection', 'VERSIONING')
    taxonomy = RelationshipFrom('Taxonomy', 'VERSIONING')
    transcriptionfactor = RelationshipFrom('TranscriptionFactor', 'VERSIONING')
    transcriptionfactorfamily = RelationshipFrom('TranscriptionFactorFamily', 'VERSIONING')
    rnafamily = RelationshipFrom('RNAFamily', 'VERSIONING')
    effector = RelationshipFrom('Effector', 'VERSIONING')
    pathway = RelationshipFrom('Pathway', 'VERSIONING')
    genome = RelationshipFrom('Genome', 'VERSIONING')
    regulog = RelationshipFrom('Regulog', 'VERSIONING')
    regulon = RelationshipFrom('Regulon', 'VERSIONING')
    operon = RelationshipFrom('Operon', 'VERSIONING')
    gene = RelationshipFrom('Gene', 'VERSIONING')

    @property
    def children_classes(self):
        return [self.details, self.collection, self.taxonomy, self.transcriptionfactor,
                self.transcriptionfactorfamily, self.rnafamily, self.effector,
                self.pathway, self.genome, self.regulog, self.regulon, self.gene]

    @property
    def children(self):

        nodes = []
        for relationship in self.children_classes:
            nodes.extend(relationship.all())

        return nodes


class BaseNode:
    uid = UniqueIdProperty()
    created = DateTimeProperty(default_now=True)

    children = ArrayProperty(IntegerProperty())

    version = RelationshipTo(Version, 'VERSIONING')

    def connect_version(self, version: Union[None, Version]):

        if version is None:
            return

        self.version.connect(version)
        rel_name = self.__class__.__name__.lower()
        node_rel = getattr(version, rel_name)
        node_rel.connect(self)


# noinspection PyAbstractClass
class Collections(StructuredNode, BaseNode):
    name = StringProperty(default='regprecise collections')
    url = StringProperty(default='https://regprecise.lbl.gov/collections.jsp')
    doi = StringProperty(default='10.1186/1471-2164-14-745')
    authors = ArrayProperty(StringProperty(), default=['Pavel S Novichkov', 'Dmitry A Rodionov'])
    description = StringProperty(default='Collection of Manually Curated Inferences of Regulons in Prokaryotic Genomes')

    children = ArrayProperty(StringProperty())


# noinspection PyAbstractClass
class Collection(StructuredNode, BaseNode):
    name = StringProperty(required=True)
    query_key = StringProperty(required=True)
    url = StringProperty()
    description = StringProperty()

    collections = RelationshipTo(Collections, 'IS_FROM')


# noinspection PyAbstractClass
class Taxonomy(StructuredNode, BaseNode):
    collection_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    genomes = ArrayProperty(IntegerProperty())

    collection = RelationshipTo(Collection, 'IS_FROM')


# noinspection PyAbstractClass
class TranscriptionFactor(StructuredNode, BaseNode):
    collection_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    full_description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    regulogs = ArrayProperty(IntegerProperty())

    collection = RelationshipTo(Collection, 'IS_FROM')


# noinspection PyAbstractClass
class TranscriptionFactorFamily(StructuredNode, BaseNode):
    tffamily_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    full_description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())

    regulogs = ArrayProperty(IntegerProperty())

    collection = RelationshipTo(Collection, 'IS_FROM')


# noinspection PyAbstractClass
class RNAFamily(StructuredNode, BaseNode):
    riboswitch_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    full_description = StringProperty()
    pubmed = ArrayProperty(IntegerProperty())
    rfam_sanger = StringProperty()

    regulogs = ArrayProperty(IntegerProperty())

    collection = RelationshipTo(Collection, 'IS_FROM')


# noinspection PyAbstractClass
class Effector(StructuredNode, BaseNode):
    effector_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    group = StringProperty()

    regulogs = ArrayProperty(IntegerProperty())

    collection = RelationshipTo(Collection, 'IS_FROM')


# noinspection PyAbstractClass
class Pathway(StructuredNode, BaseNode):
    pathway_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()
    group = StringProperty()

    regulogs = ArrayProperty(IntegerProperty())

    collection = RelationshipTo(Collection, 'IS_FROM')


# noinspection PyAbstractClass
class Genome(StructuredNode, BaseNode):
    genome_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    description = StringProperty()

    regulons = ArrayProperty(IntegerProperty())

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')


# noinspection PyAbstractClass
class Regulog(StructuredNode, BaseNode):
    regulog_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulator_type = StringProperty()
    regulator_family = StringProperty()
    regulation_mode = StringProperty()
    biological_process = StringProperty()
    regulation_effector = StringProperty()
    phylum = StringProperty()

    regulons = ArrayProperty(IntegerProperty())

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    transcription_factor = RelationshipTo(TranscriptionFactor, 'IS_FROM')
    tf_family = RelationshipTo(TranscriptionFactorFamily, 'IS_FROM')
    rna_family = RelationshipTo(RNAFamily, 'IS_FROM')
    effector = RelationshipTo(Effector, 'IS_FROM')
    pathway = RelationshipTo(Pathway, 'IS_FROM')


# noinspection PyAbstractClass
class Regulon(StructuredNode, BaseNode):
    regulon_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulator_type = StringProperty()
    locus_tag = StringProperty(required=True)
    regulator_family = StringProperty()
    regulation_mode = StringProperty()
    biological_process = StringProperty()
    regulation_effector = StringProperty()
    regulation_regulog = IntegerProperty(required=True)

    regulog = RelationshipTo(Regulog, 'IS_FROM')

    taxonomy = RelationshipTo(Taxonomy, 'IS_FROM')
    transcription_factor = RelationshipTo(TranscriptionFactor, 'IS_FROM')
    tf_family = RelationshipTo(TranscriptionFactorFamily, 'IS_FROM')
    rna_family = RelationshipTo(RNAFamily, 'IS_FROM')
    effector = RelationshipTo(Effector, 'IS_FROM')
    pathway = RelationshipTo(Pathway, 'IS_FROM')

    operons = ArrayProperty(StringProperty())
    genes = ArrayProperty(StringProperty())
    tfbs = ArrayProperty(StringProperty())


# noinspection PyAbstractClass
class Operon(StructuredNode, BaseNode):
    operon_id = IntegerProperty(required=True)
    name = StringProperty(required=True)
    url = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')
    gene = RelationshipTo('Gene', 'HAS')

    genes = ArrayProperty(StringProperty())
    tfbs = ArrayProperty(StringProperty())


# noinspection PyAbstractClass
class Gene(StructuredNode, BaseNode):
    locus_id = IntegerProperty(required=True)
    name = StringProperty()
    url = StringProperty()

    locus_tag = StringProperty(required=True)
    function = StringProperty()

    regulon = RelationshipTo(Regulon, 'IS_FROM')

    tfbs = ArrayProperty(StringProperty())


# noinspection PyAbstractClass
class TFBS(StructuredNode, BaseNode):
    regulon = RelationshipTo(Regulon, 'IS_FROM', OneOrMore)
    gene = RelationshipTo(Gene, 'IS_FROM', OneOrMore)

    position = IntegerProperty(required=True)
    score = FloatProperty(required=True)
    sequence = StringProperty(required=True)
