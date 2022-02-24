import unittest

from protrend.model import *
from protrend.utils import NeoDatabase
from protrend.utils.set_list import ProtrendSetList


class FunctionalDatabaseTest(unittest.TestCase):
    """
    Testing the database status regarding functional annotation and organisms
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

    def test_regulators_organisms_source(self):
        """
        It verifies that all regulators have one organism and sources
        """
        # noinspection PyTypeChecker
        for _regulator in Regulator.nodes:
            regulator_organism = ProtrendSetList(_regulator.organism.all())
            self.assertEqual(len(regulator_organism), 1)

            regulator_sources = ProtrendSetList(_regulator.data_source.all())
            self.assertGreaterEqual(len(regulator_sources), 1)

    def test_genes_organisms_source(self):
        """
        It verifies that all genes have one organism and sources
        """
        # noinspection PyTypeChecker
        for _gene in Gene.nodes:
            gene_organism = ProtrendSetList(_gene.organism.all())
            self.assertEqual(len(gene_organism), 1)

            gene_sources = ProtrendSetList(_gene.data_source.all())
            self.assertGreaterEqual(len(gene_sources), 1)

    def test_tfbs_organisms_source(self):
        """
        It verifies that all tfbs have one organism and sources
        """
        # noinspection PyTypeChecker
        for _tfbs in TFBS.nodes:
            tfbs_organisms = ProtrendSetList(_tfbs.data_organism.all())
            self.assertEqual(len(tfbs_organisms), 1)
            self.assertEqual(tfbs_organisms[0].protrend_id, _tfbs.organism)

            tfbs_sources = ProtrendSetList(_tfbs.data_source.all())
            self.assertGreaterEqual(len(tfbs_sources), 1)

    def test_operons_organisms_source(self):
        """
        It verifies that all operons have one organism and one source
        """
        # noinspection PyTypeChecker
        for _operon in Operon.nodes:
            operon_organisms = ProtrendSetList(_operon.organism.all())
            self.assertEqual(len(operon_organisms), 1)

            operon_sources = ProtrendSetList(_operon.data_source.all())
            self.assertEqual(len(operon_sources), 1)

    def test_regulators_locus_tag(self):
        """
        It verifies regulators with very small locus tag
        """
        # noinspection PyTypeChecker
        for _regulator in Regulator.nodes:
            self.assertGreaterEqual(len(_regulator.locus_tag), 4)

    def test_genes_locus_tag(self):
        """
        It verifies genes with very small locus tag
        """
        # noinspection PyTypeChecker
        for _gene in Gene.nodes:
            self.assertGreaterEqual(len(_gene.locus_tag), 4)

    def test_regulators_genes_organisms(self):
        """
        It verifies that regulators only regulate genes of the same organism
        """
        # noinspection PyTypeChecker
        for _regulator in Regulator.nodes:
            regulator_organisms = ProtrendSetList(_regulator.organism.all())
            self.assertEqual(len(regulator_organisms), 1)

            regulator_organism = regulator_organisms[0]
            regulator_genes = ProtrendSetList(_regulator.gene.all())

            for _gene in regulator_genes:
                gene_organism = ProtrendSetList(_gene.organism.all())[0]
                self.assertEqual(regulator_organism.protrend_id, gene_organism.protrend_id)

    def test_regulators_tfbs_organisms(self):
        """
        It verifies that regulators only bind to TFBS of the same organism
        """
        # noinspection PyTypeChecker
        for _regulator in Regulator.nodes:
            regulator_organisms = ProtrendSetList(_regulator.organism.all())
            self.assertEqual(len(regulator_organisms), 1)

            regulator_organism = ProtrendSetList(_regulator.organism.all())[0]
            regulator_tfbs = ProtrendSetList(_regulator.tfbs.all())

            for _tfbs in regulator_tfbs:
                self.assertEqual(regulator_organism.protrend_id, _tfbs.organism)

    def test_genes_tfbs_organisms(self):
        """
        It verifies that genes are only associated to TFBS of the same organism
        """
        # noinspection PyTypeChecker
        for _gene in Gene.nodes:
            gene_organism = ProtrendSetList(_gene.organism.all())[0]
            gene_tfbs = ProtrendSetList(_gene.tfbs.all())

            for _tfbs in gene_tfbs:
                self.assertEqual(gene_organism.protrend_id, _tfbs.organism)


if __name__ == '__main__':
    unittest.main()
