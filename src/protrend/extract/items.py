from itemloaders.processors import MapCompose, TakeFirst, Join
from scrapy.item import Item, Field
from w3lib.html import remove_tags

from .processors import RegPreciseProcessors


# ----------------------------
# RegPrecise items
# ----------------------------
class CollectionItem(Item):
    # outgoing relationship
    # TODO: missing add collection value
    taxonomy = Field()
    transcription_factor = Field()
    transcription_factor_family = Field()
    rna_family = Field()
    effector = Field()
    pathway = Field()


class TaxonomyItem(Item):
    collection_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                          output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    # outgoing relationship
    # TODO: missing add collection value
    collection = Field(output_processor=TakeFirst())
    genome = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class GenomeItem(Item):
    genome_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                      output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    # outgoing relationship
    taxonomy = Field(output_processor=TakeFirst())
    regulon = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class TranscriptionFactorItem(Item):
    collection_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                          output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    description = Field(input_processor=MapCompose(remove_tags),
                        output_processor=Join(separator=''))
    pubmed = Field(input_processor=MapCompose(RegPreciseProcessors.process_pubmed_href))

    # outgoing relationship
    # TODO: missing add collection value
    collection = Field(output_processor=TakeFirst())
    regulog = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class RegulogItem(Item):
    regulog_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                       output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    regulator_type = Field(output_processor=TakeFirst())

    regulator_family = Field(output_processor=TakeFirst())

    regulation_mode = Field(output_processor=TakeFirst())

    biological_process = Field(output_processor=TakeFirst())

    regulation_effector = Field(output_processor=TakeFirst())

    phylum = Field(output_processor=TakeFirst())

    # outgoing relationship
    transcription_factor = Field(output_processor=TakeFirst())
    regulon = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))

    # other outgoing relationships
    taxonomy = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    tf_family = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    rna_family = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    effector = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    pathway = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class TranscriptionFactorFamilyItem(Item):
    tffamily_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                        output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    description = Field(input_processor=MapCompose(remove_tags),
                        output_processor=Join(separator=''))
    pubmed = Field(input_processor=MapCompose(RegPreciseProcessors.process_pubmed_href))

    # outgoing relationship
    # TODO: missing add collection value
    collection = Field(output_processor=TakeFirst())
    regulog = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class RNAFamilyItem(Item):
    riboswitch_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                          output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    description = Field(input_processor=MapCompose(remove_tags),
                        output_processor=Join(separator=''))
    pubmed = Field(input_processor=MapCompose(RegPreciseProcessors.process_pubmed_href))
    rfam = Field(input_processor=MapCompose(RegPreciseProcessors.process_rfam_href),
                 output_processor=TakeFirst())

    # outgoing relationship
    # TODO: missing add collection value
    collection = Field(output_processor=TakeFirst())
    regulog = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class EffectorItem(Item):
    effector_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                        output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    # outgoing relationship
    # TODO: missing add collection value
    collection = Field(output_processor=TakeFirst())
    regulog = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class PathwayItem(Item):
    pathway_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                       output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    # outgoing relationship
    # TODO: missing add collection value
    collection = Field(output_processor=TakeFirst())
    regulog = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class RegulonItem(Item):
    regulon_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_href),
                       output_processor=TakeFirst())

    name = Field(output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    regulator_type = Field(output_processor=TakeFirst())

    regulator_locus_tag = Field(output_processor=TakeFirst())

    regulator_family = Field(output_processor=TakeFirst())

    regulation_mode = Field(output_processor=TakeFirst())

    biological_process = Field(output_processor=TakeFirst())

    regulation_effector = Field(output_processor=TakeFirst())

    regulation_regulog = Field(input_processor=MapCompose(RegPreciseProcessors.process_regulog_name),
                               output_processor=TakeFirst())

    # outgoing relationship
    genome = Field(output_processor=TakeFirst())
    operon = Field()
    gene = Field()
    tfbs = Field()

    # other outgoing relationships
    regulog = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    taxonomy = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    transcription_factor = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    tf_family = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    rna_family = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    effector = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))
    pathway = Field(input_processor=MapCompose(RegPreciseProcessors.process_href))


class OperonItem(Item):
    operon_id = Field(output_processor=Join(separator='_'))

    name = Field(input_processor=MapCompose(RegPreciseProcessors.process_operon_name),
                 output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    # outgoing relationship
    regulon = Field(output_processor=TakeFirst())
    gene = Field()
    tfbs = Field()


class GeneItem(Item):
    locus_tag = Field(input_processor=MapCompose(RegPreciseProcessors.process_locus_tag),
                      output_processor=TakeFirst())

    name = Field(input_processor=MapCompose(RegPreciseProcessors.process_name),
                 output_processor=TakeFirst())

    function = Field(input_processor=MapCompose(RegPreciseProcessors.process_function),
                     output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    # outgoing relationship
    regulon = Field(output_processor=TakeFirst())
    operon = Field(output_processor=TakeFirst())
    tfbs = Field()


class TFBSItem(Item):
    tfbs_id = Field(input_processor=MapCompose(RegPreciseProcessors.process_score_str),
                    output_processor=Join(separator='_'))

    position = Field(input_processor=MapCompose(RegPreciseProcessors.process_position),
                     output_processor=TakeFirst())

    score = Field(input_processor=MapCompose(RegPreciseProcessors.process_score),
                  output_processor=TakeFirst())

    sequence = Field(input_processor=MapCompose(RegPreciseProcessors.process_sequence),
                     output_processor=TakeFirst())

    url = Field(output_processor=TakeFirst())

    # outgoing relationship
    regulon = Field(output_processor=TakeFirst())
    operon = Field(output_processor=TakeFirst())
    gene = Field()
