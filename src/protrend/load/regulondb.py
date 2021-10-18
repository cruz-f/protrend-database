from protrend.load.loader import Loader


class RegulondbLoader(Loader):
    default_source = 'regulondb'
    default_version = '0.0.0'
    default_load_stack = [
        'nodes_effector.json',
        'nodes_evidence.json',
        'nodes_gene.json',
        'nodes_operon.json',
        'nodes_organism.json',
        'nodes_promoter.json',
        'nodes_publication.json',
        'nodes_regulator.json',
        'nodes_regulatoryfamily.json',
        'nodes_regulatoryinteraction.json',
        'nodes_source.json',
        'nodes_tfbs.json',

        'connected_evidence_gene.json',
        'connected_evidence_operon.json',
        'connected_evidence_promoter.json',
        'connected_evidence_regulator.json',
        'connected_evidence_regulatoryinteraction.json',
        'connected_evidence_tfbs.json',

        'connected_gene_promoter.json',
        'connected_gene_tfbs.json',
        'connected_operon_gene.json',
        'connected_operon_promoter.json',
        'connected_operon_tfbs.json',

        'connected_regulator_gene.json',
        'connected_regulator_operon.json',
        'connected_regulator_tfbs.json',

        'connected_effector_organism.json',
        'connected_gene_organism.json',
        'connected_operon_organism.json',
        'connected_promoter_organism.json',
        'connected_regulator_organism.json',
        'connected_regulatoryinteraction_organism.json',
        'connected_tfbs_organism.json',

        'connected_publication_gene.json',
        'connected_publication_operon.json',
        'connected_publication_promoter.json',
        'connected_publication_regulator.json',
        'connected_publication_tfbs.json',

        'connected_regulator_regulatoryfamily.json',

        'connected_regulator_effector.json',
        'connected_regulator_gene.json',
        'connected_regulator_operon.json',
        'connected_regulator_tfbs.json',

        'connected_regulatoryinteraction_effector.json',
        'connected_regulatoryinteraction_gene.json',
        'connected_regulatoryinteraction_operon.json',
        'connected_regulatoryinteraction_regulator.json',
        'connected_regulatoryinteraction_tfbs.json',

        'connected_effector_source.json',
        'connected_gene_source.json',
        'connected_operon_source.json',
        'connected_organism_source.json',
        'connected_promoter_source.json',
        'connected_regulator_source.json',
        'connected_regulatoryfamily_source.json',
        'connected_regulatoryinteraction_source.json',
        'connected_tfbs_source.json',
    ]
