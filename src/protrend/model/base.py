from collections import defaultdict
from datetime import datetime
from typing import List, Dict, Any, Union, Type

import pandas as pd
import pytz
from neomodel import (UniqueIdProperty, DateTimeProperty, StructuredNode, StringProperty, RelationshipManager,
                      IntegerProperty, ArrayProperty)
from tqdm import tqdm

from protrend.utils.miscellaneous import convert_to_snake_case, is_null
from .utils import choices, help_text, _sort_nodes


class BaseNode(StructuredNode):
    __abstract_node__ = True

    uid = UniqueIdProperty()
    protrend_id = StringProperty(required=True, unique_index=True, help_text=help_text.protrend_id)
    created = DateTimeProperty(default_now=True, help_text=help_text.created)
    updated = DateTimeProperty(default_now=True, help_text=help_text.updated)

    identifying_property = 'protrend_id'
    header = 'PRT'
    entity = 'PRT'
    node_factors = {}

    node_register: Dict[str, Type['BaseNode']] = {}

    def __init_subclass__(cls, **kwargs):
        cls.node_register[cls.node_name()] = cls

    # -------------------------------------
    # Class attributes
    # -------------------------------------
    @classmethod
    def node_name(cls):
        return convert_to_snake_case(cls.__name__.lower())

    @classmethod
    def node_properties(cls) -> dict:
        return dict(cls.__all_properties__)

    @classmethod
    def node_relationships(cls) -> Dict[str, RelationshipManager]:
        return dict(cls.__all_relationships__)

    @classmethod
    def node_keys(cls):
        return cls.node_properties().keys()

    @classmethod
    def node_values(cls):
        return cls.node_properties().values()

    @classmethod
    def node_items(cls):
        return cls.node_properties().items()

    @classmethod
    def last_node(cls) -> Union['BaseNode', None]:
        nodes = cls.nodes.all()
        if nodes:
            sorted_nodes = sorted(nodes, key=_sort_nodes, reverse=True)
            return sorted_nodes[0]

        return None

    # -------------------------------------
    # Class methods
    # -------------------------------------
    @classmethod
    def node_from_dict(cls,
                       *nodes: Dict[str, Any],
                       save: bool = False) -> List['BaseNode']:

        structured_nodes = []
        node_keys = list(cls.node_keys())

        for node in nodes:

            node_kwargs = {key: val for key, val in node.items()
                           if key in node_keys and val is not None}

            if node_kwargs:

                structured_node = cls(**node_kwargs)

                if save:
                    structured_node = structured_node.save()

                structured_nodes.append(structured_node)

        return structured_nodes

    @classmethod
    def node_update_from_dict(cls,
                              *nodes: Dict[str, Any],
                              save: bool = False) -> List['BaseNode']:

        structured_nodes = []
        node_keys = [key for key in cls.node_keys() if key != cls.identifying_property]
        all_nodes = cls.node_to_dict(to='node')

        for node in nodes:

            identifier = node.get(cls.identifying_property, '')

            if identifier in all_nodes:

                structured_node = all_nodes[identifier]

                for key, val in node.items():

                    if key in node_keys:

                        last_val = getattr(structured_node, key, None)

                        if not is_null(val) and is_null(last_val):
                            setattr(structured_node, key, val)

                if save:
                    structured_node.updated = datetime.utcnow().replace(tzinfo=pytz.utc)
                    structured_node = structured_node.save()

                structured_nodes.append(structured_node)

        return structured_nodes

    @classmethod
    def node_to_dict(cls, to: str = 'dict') -> Union[Dict[str, Any], Dict[str, 'BaseNode']]:

        if to == 'dict':
            res = defaultdict(list)

            for node in cls.nodes.all():

                node_properties = node.properties

                for key in cls.node_keys():
                    val = node_properties.get(key, None)
                    res[key].append(val)

            if not res:
                return {key: [] for key in cls.node_keys()}

            return res

        elif to == 'node':

            return {node.identifier: node for node in cls.nodes.all()}

        raise ValueError(f'Invalid output {to}')

    @classmethod
    def node_from_df(cls,
                     nodes: pd.DataFrame,
                     save: bool = False) -> List['BaseNode']:

        structured_nodes = []
        node_keys = list(cls.node_keys())
        node_name = cls.node_name()

        for _, node in tqdm(nodes.iterrows(), desc=f'{node_name} - node_creation', total=nodes.shape[0]):

            node_kwargs = {key: val for key, val in node.items()
                           if key in node_keys and not is_null(val)}

            if node_kwargs:

                structured_node = cls(**node_kwargs)

                if save:
                    structured_node = structured_node.save()

                structured_nodes.append(structured_node)

        return structured_nodes

    @classmethod
    def node_update_from_df(cls,
                            nodes: pd.DataFrame,
                            save: bool = False) -> List['BaseNode']:

        structured_nodes = []
        node_keys = [key for key in cls.node_keys() if key != cls.identifying_property]
        all_nodes = cls.node_to_dict(to='node')
        node_name = cls.node_name()

        for _, node in tqdm(nodes.iterrows(), desc=f'{node_name} - node_update', total=nodes.shape[0]):

            identifier = node.get(cls.identifying_property, '')

            if identifier in all_nodes:

                structured_node = all_nodes[identifier]

                for key, val in node.items():

                    key: str

                    if key in node_keys:

                        last_val = getattr(structured_node, key, None)

                        if not is_null(val) and is_null(last_val):
                            setattr(structured_node, key, val)

                if save:
                    structured_node.updated = datetime.utcnow().replace(tzinfo=pytz.utc)
                    structured_node = structured_node.save()

                structured_nodes.append(structured_node)

        return structured_nodes

    @classmethod
    def node_to_df(cls) -> pd.DataFrame:
        return pd.DataFrame.from_dict(cls.node_to_dict())

    # -------------------------------------
    # Instance methods
    # -------------------------------------
    @property
    def identifier(self):
        return self.properties[self.identifying_property]

    @property
    def properties(self) -> Dict[str, Any]:
        return {key: val for key, val in self.__properties__.items()
                if val is not None and key != 'id'}

    def keys(self):
        return self.properties.keys()

    def values(self):
        return self.properties.values()

    def items(self):
        return self.properties.items()

    def to_dict(self):
        return self.properties

    def to_series(self):
        return pd.Series(data=self.values(), index=self.keys())

    def to_df(self):
        data = [list(self.values())]
        return pd.DataFrame(data=data, columns=self.keys())


class NameMixIn:
    # properties
    name = StringProperty(required=True, unique_index=True, max_length=250, help_text=help_text.required_name)
    name_factor = StringProperty(required=True, unique_index=True, max_length=250, help_text=help_text.required_name)


class SequenceMixIn:
    sequence = StringProperty(help_text=help_text.sequence)


class PositionMixIn:
    strand = StringProperty(choices=choices.strand, help_text=help_text.strand)
    start = IntegerProperty(help_text=help_text.start)
    stop = IntegerProperty(help_text=help_text.stop)


class GeneMixIn:
    # properties
    locus_tag = StringProperty(required=True, unique_index=True, max_length=50, help_text=help_text.locus_tag)
    locus_tag_factor = StringProperty(required=True, unique_index=True, max_length=50,
                                      help_text=help_text.required_name)
    uniprot_accession = StringProperty(max_length=50, help_text=help_text.uniprot_accession)
    uniprot_accession_factor = StringProperty(max_length=50, help_text=help_text.uniprot_accession)
    name = StringProperty(max_length=50, help_text=help_text.gene_name)
    synonyms = ArrayProperty(StringProperty(), help_text=help_text.synonyms)
    function = StringProperty(help_text=help_text.function)
    description = StringProperty(help_text=help_text.description)
    ncbi_gene = IntegerProperty(max_length=50, help_text=help_text.ncbi_gene)
    ncbi_protein = IntegerProperty(max_length=50, help_text=help_text.ncbi_protein)
    genbank_accession = StringProperty(max_length=50, help_text=help_text.genbank_accession)
    refseq_accession = StringProperty(max_length=50, help_text=help_text.refseq_accession)
