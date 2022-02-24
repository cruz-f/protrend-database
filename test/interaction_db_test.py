import unittest

from protrend.model import *
from protrend.utils import NeoDatabase
from protrend.utils.set_list import ProtrendSetList


class InteractionDatabaseTest(unittest.TestCase):
    """
    Testing the database status regarding functional annotation and interactions
    """

    def setUp(self) -> None:
        """
        Set up to connect to the database
        :return:
        """
        user_name = 'neo4j'
        password = 'protrend'
        ip = 'localhost'
        port = '7687'
        neo_db = NeoDatabase(user_name=user_name, password=password, ip=ip, port=port)
        neo_db.connect()

    def test_interactions(self):
        """
        It verifies that all interactions have the correct children
        """
        # noinspection PyTypeChecker
        for interaction in RegulatoryInteraction.nodes:
            interaction_sources = ProtrendSetList(interaction.data_source.all())
            self.assertGreaterEqual(len(interaction_sources), 1)

            interaction_organisms = ProtrendSetList(interaction.data_organism.all())
            interaction_regulators = ProtrendSetList(interaction.data_regulator.all())
            interaction_genes = ProtrendSetList(interaction.data_gene.all())
            interaction_tfbs = ProtrendSetList(interaction.data_tfbs.all())
            interaction_effectors = ProtrendSetList(interaction.data_effector.all())

            # interaction - organisms
            self.assertEqual(len(interaction_organisms), 1)
            self.assertEqual(interaction_organisms[0].protrend_id, interaction.organism)

            # interaction - regulators
            self.assertEqual(len(interaction_regulators), 1)
            self.assertEqual(interaction_regulators[0].protrend_id, interaction.regulator)

            # interaction - genes
            self.assertEqual(len(interaction_genes), 1)
            self.assertEqual(interaction_genes[0].protrend_id, interaction.gene)

            # interaction - tfbs
            self.assertEqual(len(interaction_tfbs), 1)
            self.assertEqual(interaction_tfbs[0].protrend_id, interaction.tfbs)

            # interaction - effectors
            self.assertEqual(len(interaction_effectors), 1)
            self.assertEqual(interaction_effectors[0].protrend_id, interaction.effector)

    def test_interactions_regulators(self):
        """
        It verifies that all interactions have the correct children
        """
        # noinspection PyTypeChecker
        for interaction in RegulatoryInteraction.nodes:

            interaction_regulators = ProtrendSetList(interaction.data_regulator.all())
            _regulator = interaction_regulators[0]

            # regulator - organisms
            regulator_organism = _regulator.organism[0]
            self.assertEqual(interaction.organism, regulator_organism.protrend_id)

            # regulator - genes
            regulator_genes = [obj.protrend_id for obj in _regulator.gene.all()]
            self.assertIn(interaction.gene, regulator_genes)

            # regulator - tfbs
            regulator_tfbs = [obj.protrend_id for obj in _regulator.tfbs.all()]
            self.assertIn(interaction.tfbs, regulator_tfbs)

            # regulator - effectors
            regulator_effectors = [obj.protrend_id for obj in _regulator.effector.all()]
            self.assertIn(interaction.effector, regulator_effectors)

    def test_interactions_genes(self):
        """
        It verifies that all interactions have the correct children
        """
        # noinspection PyTypeChecker
        for interaction in RegulatoryInteraction.nodes:
            interaction_genes = ProtrendSetList(interaction.data_gene.all())
            _gene = interaction_genes[0]

            # gene - organisms
            gene_organism = _gene.organism[0]
            self.assertEqual(interaction.organism, gene_organism.protrend_id)

            # gene - regulators
            gene_regulators = [obj.protrend_id for obj in _gene.regulator.all()]
            self.assertIn(interaction.regulator, gene_regulators)

            # gene - tfbs
            gene_tfbs = [obj.protrend_id for obj in _gene.tfbs.all()]
            self.assertIn(interaction.tfbs, gene_tfbs)

    def test_interactions_tfbs(self):
        """
        It verifies that all interactions have the correct children
        """
        # noinspection PyTypeChecker
        for interaction in RegulatoryInteraction.nodes:
            interaction_tfbs = ProtrendSetList(interaction.data_tfbs.all())
            _tfbs = interaction_tfbs[0]

            # tfbs - organisms
            tfbs_organism = _tfbs.data_organism[0]
            self.assertEqual(interaction.organism, tfbs_organism.protrend_id)

            # tfbs - regulators
            tfbs_regulators = [obj.protrend_id for obj in _tfbs.regulator.all()]
            self.assertIn(interaction.regulator, tfbs_regulators)

            # tfbs - genes
            tfbs_genes = [obj.protrend_id for obj in _tfbs.gene.all()]
            self.assertIn(interaction.tfbs, tfbs_genes)

    def test_interactions_effectors(self):
        """
        It verifies that all interactions have the correct children
        """
        # noinspection PyTypeChecker
        for interaction in RegulatoryInteraction.nodes:
            interaction_effector = ProtrendSetList(interaction.data_effector.all())
            _effector = interaction_effector[0]

            # effector - regulators
            effector_regulators = [obj.protrend_id for obj in _effector.regulator.all()]
            self.assertIn(interaction.regulator, effector_regulators)

    def test_interactions_regulators_genes(self):
        """
        It verifies that all interactions have the correct children
        """
        # noinspection PyTypeChecker
        for _regulator in Regulator.nodes:
            interactions = RegulatoryInteraction.nodes.filter(regulator__exact=_regulator.protrend_id)

            interactions_genes = {obj.gene for obj in interactions}
            regulator_genes = {obj.protrend_id for obj in _regulator.gene.all()}

            self.assertEqual(len(interactions_genes), len(regulator_genes))

            for interaction_gene in interactions_genes:
                self.assertIn(interaction_gene, regulator_genes)


if __name__ == '__main__':
    unittest.main()
