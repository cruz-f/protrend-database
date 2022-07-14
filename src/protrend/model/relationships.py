from neomodel import StructuredRel, DateTimeProperty, StringProperty, IntegerProperty

from protrend.model.utils import help_text
from protrend.model.utils import choices


BASE_REL_TYPE = 'HAS'
SOURCE_REL_TYPE = 'OWNER'
ALIGNED_SEQUENCE_REL_TYPE = 'ALIGNED_SEQUENCE'


class BaseRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)


class SourceRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    key = StringProperty()
    url = StringProperty()
    external_identifier = StringProperty()


class AlignedSequenceRelationship(StructuredRel):
    # base
    created = DateTimeProperty(default_now=True)
    updated = DateTimeProperty(default_now=True)

    # properties
    sequence = StringProperty(help_text=help_text.aligned_sequence)
    start = IntegerProperty(help_text=help_text.start)
    stop = IntegerProperty(help_text=help_text.stop)
    strand = StringProperty(choices=choices.strand, help_text=help_text.strand)
