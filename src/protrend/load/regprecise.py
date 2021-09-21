from protrend.load.loader import Loader


class RegPreciseLoader(Loader):
    default_source = 'regprecise'
    default_version = '0.0.0'
    default_load_stack = [
        'nodes_effector.json',
        'nodes_gene.json',
        'nodes_operon.json',
        'nodes_organism.json',
        'nodes_pathway.json',
        'nodes_publication.json',
        'nodes_regulator.json',
        'nodes_regulatoryfamily.json',
        'nodes_regulatoryinteraction.json',
        'nodes_source.json',
        'nodes_tfbs.json',

        'connected_effector_organism.json',
        'connected_effector_regulator.json',
        'connected_effector_source.json',

        'connected_gene_source.json',
        'connected_gene_organism.json',

        'connected_gene_regulator.json',
        'connected_gene_tfbs.json',

        'connected_operon_gene.json',
        'connected_operon_organism.json',
        'connected_operon_regulator.json',
        'connected_operon_source.json',
        'connected_operon_tfbs.json',
        'connected_tfbs_regulator.json',

        'connected_organism_source.json',

        'connected_pathway_gene.json',
        'connected_pathway_regulator.json',
        'connected_pathway_source.json',

        'connected_regulator_organism.json',
        'connected_regulator_source.json',

        'connected_regulatoryfamily_publication.json',
        'connected_regulatoryfamily_regulator.json',
        'connected_regulatoryfamily_source.json',

        'connected_regulatoryinteraction_effector.json',
        'connected_regulatoryinteraction_gene.json',
        'connected_regulatoryinteraction_operon.json',
        'connected_regulatoryinteraction_organism.json',
        'connected_regulatoryinteraction_regulator.json',
        'connected_regulatoryinteraction_source.json',
        'connected_regulatoryinteraction_tfbs.json',

        'connected_tfbs_organism.json',
        'connected_tfbs_source.json',
    ]
