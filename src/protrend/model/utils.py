from typing import Union, Type, Tuple, List

from neomodel import RelationshipManager, StructuredRel, StringProperty

from protrend import BaseNode

FORWARD = 'FORWARD'
REVERSE = 'REVERSE'
UNKNOWN = 'UNKNOWN'
TRANSCRIPTION_FACTOR = 'TRANSCRIPTION FACTOR'
TRANSCRIPTION_ATTENUATOR = 'TRANSCRIPTION ATTENUATOR'
TRANSCRIPTION_TERMINATOR = 'TRANSCRIPTION TERMINATOR'
SIGMA_FACTOR = 'SIGMA FACTOR'
SMALL_RNA = 'SMALL RNA (sRNA)'
ACTIVATION = 'ACTIVATION'
REPRESSION = 'REPRESSION'
DUAL = 'DUAL'
LITERATURE = 'LITERATURE'
DATABASE = 'DATABASE'
CURATION = 'CURATION'


# noinspection PyPep8Naming
class help_text:
    protrend_id = 'Universal identifier of the ProTReND database'
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
    sequence = 'The protein sequence for this protein'
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
    effector_id = 'The effector ProTReND identifier'
    regulatory_effect = 'The regulatory effect (eg activation, repression, dual and unknown) of this regulatory interaction'
    source_author = 'The authors of this data source'
    url = 'The web address for this data source'
    data_source_type = 'The type of this data source'
    length = 'The length of the TFBS sequence'


# noinspection PyPep8Naming
class choices:
    strand = {FORWARD: 'forward', REVERSE: 'reverse', UNKNOWN: 'unknown'}
    mechanism = {TRANSCRIPTION_FACTOR: 'transcription factor', TRANSCRIPTION_ATTENUATOR: 'transcription attenuator',
                 TRANSCRIPTION_TERMINATOR: 'transcription terminator', SIGMA_FACTOR: 'sigma factor',
                 SMALL_RNA: 'small RNA (sRNA)', UNKNOWN: 'unknown'}
    regulatory_effect = {ACTIVATION: 'activation', REPRESSION: 'repression', DUAL: 'dual', UNKNOWN: 'unknown'}
    data_source_type = {LITERATURE: 'literature', DATABASE: 'database', CURATION: 'curation'}


def _sort_nodes(node: BaseNode):
    return protrend_id_decoder(node.protrend_id)


def protrend_id_encoder(header: str, entity: str, integer: Union[str, int]) -> str:
    integer = int(integer)

    return f'{header}.{entity}.{integer:07}'


def protrend_id_decoder(protrend_id: Union[str, StringProperty]) -> int:
    prt, entity, integer = protrend_id.split('.')

    return int(integer)


def get_node_by_name(name: str, default=None) -> Union[Type[BaseNode], None]:
    return BaseNode.node_register.get(name, default)


def _find_to_node(relationship: RelationshipManager) -> Type[BaseNode]:
    return relationship.definition['node_class']


def get_nodes_relationships(from_node: Type[BaseNode], to_node: Type[BaseNode]) -> Tuple[List[str], List[str]]:
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


def connect_nodes(from_node: BaseNode, to_node: BaseNode, relationship: str, kwargs: dict) -> bool:
    relationship: RelationshipManager = getattr(from_node, relationship, None)
    relationship_model: StructuredRel = relationship.definition['model']

    if kwargs:
        kwargs = {key: val for key, val in kwargs.items()
                  if hasattr(relationship_model, key)}

    else:
        kwargs = {}

    return relationship.connect(to_node, kwargs)
