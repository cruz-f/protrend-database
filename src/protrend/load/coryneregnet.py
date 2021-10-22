from protrend.load.loader import Loader


class CoryneRegNetLoader(Loader):
    default_source = 'coryneregnet'
    default_version = '0.0.0'
    default_load_stack = [

        'nodes_evidence.json',
        'nodes_gene.json',
        'nodes_operon.json',
        'nodes_organism.json',
        'nodes_publication.json',
        'nodes_regulator.json',
        'nodes_regulatoryinteraction.json',
        'nodes_source.json',
        'nodes_tfbs.json',

        'connected_evidence_regulatoryinteraction.json',

        'connected_gene_tfbs.json',
        'connected_operon_gene.json',
        'connected_operon_tfbs.json',

        'connected_gene_organism.json',
        'connected_operon_organism.json',
        'connected_regulator_organism.json',
        'connected_regulatoryinteraction_organism.json',
        'connected_tfbs_organism.json',

        'connected_publication_gene.json',
        'connected_publication_operon.json',
        'connected_publication_organism.json',
        'connected_publication_regulator.json',
        'connected_publication_tfbs.json',
        'connected_publication_regulatoryinteraction.json',

        'connected_regulator_gene.json',
        'connected_regulator_operon.json',
        'connected_regulator_tfbs.json',

        'connected_regulatoryinteraction_gene.json',
        'connected_regulatoryinteraction_operon.json',
        'connected_regulatoryinteraction_regulator.json',
        'connected_regulatoryinteraction_tfbs.json',

        'connected_gene_source.json',
        'connected_operon_source.json',
        'connected_organism_source.json',
        'connected_regulator_source.json',
        'connected_regulatoryinteraction_source.json',
        'connected_tfbs_source.json',
    ]
