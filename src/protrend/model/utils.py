from typing import Union, Type, Tuple, List, TYPE_CHECKING

from neomodel import RelationshipManager, StructuredRel, StringProperty

import protrend.utils.constants as constants

if TYPE_CHECKING:
    from .base import BaseNode


# noinspection PyPep8Naming
class help_text:
    protrend_id = 'Universal identifier of the ProTReND database'
    created = 'Time tag for the item creation'
    updated = 'Time tag of the item last alteration'
    synonyms = 'A list of synonyms for this gene'
    locus_tag = 'The locus tag for this gene'
    uniprot_accession = 'The UniProt accession for this protein'
    gene_name = 'The name for this gene/protein'
    function = 'The function for this protein'
    description = 'The description for this protein'
    ncbi_gene = 'The NCBI gene identifier'
    ncbi_protein = 'The NCBI protein identifier'
    genbank_accession = 'The NCBI GenBank accession'
    refseq_accession = 'The NCBI RefSeq accession'
    protein_sequence = 'The protein sequence for this protein'
    gene_sequence = 'The gene sequence for this gene'
    strand = 'The strand corresponds to the genomic coordinate forward or reverse'
    start = 'The start corresponds to the genomic coordinate of the item position in the genome sequence'
    stop = 'The stop corresponds to the genomic coordinate of the item position in the genome sequence'
    required_name = 'The name for this item which will be used as main identifier'
    mechanism = 'The regulatory mechanism associated to this regulator'
    organism_name = 'The scientific name for this organism including strain name whenever possible'
    ncbi_taxonomy = 'The NCBI taxonomy identifier'
    species = 'The scientific name for this species'
    strain = 'The strain for this species'
    refseq_ftp = 'The NCBI RefSeq ftp address for this accession'
    genbank_ftp = 'The NCBI GenBank ftp address for this accession'
    ncbi_assembly = 'The NCBI Assembly identifier for the organism genome sequence'
    assembly_accession = 'The NCBI Assembly accession for the organism genome sequence'
    kegg_compounds = 'A list of KEGG compound identifiers associated with this effector'
    generic_description = 'The description for this item'
    operon_db_id = 'The OperonDB identifier for this operon'
    operon_name = 'The name for this operon'
    operon_function = 'The function for this operon'
    operon_genes = 'The identifiers for the genes associated with this operon'
    kegg_pathways = 'A list of KEGG pathway identifiers associated with this pathway'
    pmid = 'The PubMed identifier for this publication'
    doi = 'The Digital Object Identifier (DOI) for this publication'
    title = 'The title of this publication'
    author = 'The main author of this publication'
    year = 'The year of this publication'
    rfam = 'The Regulatory Family (RFAM) name'
    organism_id = 'The organism ProTReND identifier'
    regulator_id = 'The regulator ProTReND identifier'
    gene_id = 'The gene ProTReND identifier'
    tfbs_id = 'The TFBS ProTReND identifier'
    tfbs_sequence = 'The binding sequence for this TFBS'
    effector_id = 'The effector ProTReND identifier'
    regulatory_effect = 'The regulatory effect (eg activation, repression, dual and unknown) of this regulatory interaction'
    source_author = 'The authors of this data source'
    url = 'The web address for this data source'
    data_source_type = 'The type of this data source'
    length = 'The length of the TFBS sequence'
    consensus_sequence = 'The consensus sequence of the motif'
    aligned_sequences = 'The aligned sequences of the motif'
    aligned_sequence = 'An aligned sequence TFBS'
    promoter_sequence = 'The sequence of the promoter'


# noinspection PyPep8Naming
class choices:
    strand = {constants.FORWARD: 'forward', constants.REVERSE: 'reverse', constants.UNKNOWN: 'unknown'}
    mechanism = {constants.TRANSCRIPTION_FACTOR: 'transcription factor',
                 constants.TRANSCRIPTION_ATTENUATOR: 'transcription attenuator',
                 constants.TRANSCRIPTION_TERMINATOR: 'transcription terminator',
                 constants.SIGMA_FACTOR: 'sigma factor',
                 constants.SMALL_RNA: 'small RNA (sRNA)',
                 constants.UNKNOWN: 'unknown'}
    regulatory_effect = {constants.ACTIVATION: 'activation',
                         constants.REPRESSION: 'repression',
                         constants.DUAL: 'dual',
                         constants.UNKNOWN: 'unknown'}
    data_source_type = {constants.LITERATURE: 'literature',
                        constants.DATABASE: 'database',
                        constants.CURATION: 'curation'}


def _sort_nodes(node: 'BaseNode'):
    return protrend_id_decoder(node.protrend_id)


def protrend_id_encoder(header: str, entity: str, integer: Union[str, int]) -> str:
    integer = int(integer)

    return f'{header}.{entity}.{integer:07}'


def protrend_id_decoder(protrend_id: Union[str, StringProperty]) -> int:
    prt, entity, integer = protrend_id.split('.')

    return int(integer)


def get_node_by_name(name: str, default=None) -> Union[Type['BaseNode'], None]:
    from protrend.model.base import BaseNode
    return BaseNode.node_register.get(name, default)


def _find_to_node(relationship) -> Type['BaseNode']:
    if 'node_class' not in relationship.definition:
        # noinspection PyProtectedMember
        relationship._lookup_node_class()
    return relationship.definition['node_class']


def get_nodes_relationships(from_node: Type['BaseNode'], to_node: Type['BaseNode']) -> Tuple[List[str], List[str]]:
    from_node_rels = from_node.node_relationships()
    from_node_matches = []
    for attr, relationship in from_node_rels.items():

        this_to_node = _find_to_node(relationship)
        if this_to_node.node_name() == to_node.node_name():
            from_node_matches.append(attr)

    to_node_rels = to_node.node_relationships()
    to_node_matches = []
    for attr, relationship in to_node_rels.items():

        this_from_node = _find_to_node(relationship)
        if this_from_node.node_name() == from_node.node_name():
            to_node_matches.append(attr)

    return from_node_matches, to_node_matches


def connect_nodes(from_node: 'BaseNode',
                  to_node: 'BaseNode',
                  relationship: str,
                  kwargs: dict,
                  avoid_duplicates: bool = False) -> bool:
    relationship: RelationshipManager = getattr(from_node, relationship, None)
    relationship_model: StructuredRel = relationship.definition['model']

    if kwargs:
        kwargs = {key: val for key, val in kwargs.items()
                  if hasattr(relationship_model, key)}

    else:
        kwargs = {}

    if avoid_duplicates and relationship.is_connected(to_node):
        return False

    return relationship.connect(to_node, kwargs)
