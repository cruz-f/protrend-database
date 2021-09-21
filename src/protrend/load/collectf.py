from protrend.load.loader import Loader


class CollectfLoader(Loader):
    default_source = 'collectf'
    default_version = '0.0.1'
    default_load_stack = [
        'nodes_evidence.json',
        'nodes_gene.json',
        'nodes_operon.json',
        'nodes_organism.json',
        'nodes_publication.json',
        'nodes_regulator.json',
        'nodes_regulatoryfamily.json',
        'nodes_regulatoryinteraction.json',
        'nodes_source.json',
        'nodes_tfbs.json',

        'connected_evidence_gene.json',
        'connected_evidence_operon.json',
        'connected_evidence_regulator.json',
        'connected_evidence_regulatoryinteraction.json',
        'connected_evidence_tfbs.json',

        'connected_gene_tfbs.json',
        'connected_operon_gene.json',
        'connected_operon_tfbs.json',
        'connected_regulator_gene.json',
        'connected_regulator_operon.json',
        'connected_regulator_tfbs.json',

        'connected_organism_gene.json',
        'connected_organism_operon.json',
        'connected_organism_regulator.json',
        'connected_organism_regulatoryinteraction.json',
        'connected_organism_tfbs.json',

        'connected_publication_gene.json',
        'connected_publication_operon.json',
        'connected_publication_regulator.json',
        'connected_publication_regulatoryinteraction.json',
        'connected_publication_tfbs.json',

        'connected_regulatoryfamily_regulator.json',

        'connected_regulatoryinteraction_gene.json',
        'connected_regulatoryinteraction_operon.json',
        'connected_regulatoryinteraction_regulator.json',
        'connected_regulatoryinteraction_tfbs.json',

        'connected_gene_source.json',
        'connected_operon_source.json',
        'connected_organism_source.json',
        'connected_regulator_source.json',

        'connected_regulatoryfamily_source.json',
        'connected_regulatoryinteraction_source.json',
        'connected_tfbs_source.json',

    ]