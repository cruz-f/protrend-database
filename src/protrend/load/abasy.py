from protrend.load.loader import Loader


class AbasyLoader(Loader):
    default_source = 'abasy'
    default_version = '0.0.0'
    default_load_stack = [

        'nodes_gene.json',
        'nodes_operon.json',
        'nodes_organism.json',
        'nodes_regulator.json',
        'nodes_regulatoryinteraction.json',
        'nodes_source.json',

        'connected_operon_gene.json',

        'connected_gene_organism.json',
        'connected_operon_organism.json',
        'connected_regulator_organism.json',
        'connected_regulatoryinteraction_organism.json',

        'connected_regulator_gene.json',
        'connected_regulator_operon.json',

        'connected_regulatoryinteraction_gene.json',
        'connected_regulatoryinteraction_operon.json',
        'connected_regulatoryinteraction_regulator.json',

        'connected_gene_source.json',
        'connected_operon_source.json',
        'connected_organism_source.json',
        'connected_regulator_source.json',
        'connected_regulatoryinteraction_source.json',
    ]
