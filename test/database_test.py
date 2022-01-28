import unittest

from protrend.model import *
from protrend.utils import NeoDatabase


class DatabaseTest(unittest.TestCase):
    """
    Testing the database status
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

    def test_effector(self):
        """
        Test Effector nodes
        :return:
        """
        df = Effector.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 500)
        self.assertEqual(properties, 7)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.name.nunique(), entries)

    def test_evidence(self):
        """
        Test Evidence nodes
        :return:
        """
        df = Evidence.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 180)
        self.assertEqual(properties, 7)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.name.nunique(), entries)

    def test_gene(self):
        """
        Test Gene nodes
        :return:
        """
        df = Gene.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 114000)
        self.assertEqual(properties, 20)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.locus_tag.nunique(), entries)

        uniprot_accessions = df.uniprot_accession.nunique()
        self.assertNotEqual(uniprot_accessions, entries)
        self.assertEqual(uniprot_accessions, df.uniprot_accession.count())
        self.assertGreaterEqual(uniprot_accessions / entries, 0.75)

    def test_operon(self):
        """
        Test Operon nodes
        :return:
        """
        df = Operon.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 26000)
        self.assertEqual(properties, 12)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.operon_db_id.nunique(), entries)

        operon_by_gene = df.explode('genes')
        self.assertNotEqual(operon_by_gene.protrend_id.nunique(), operon_by_gene.genes.nunique())

    def test_organism(self):
        """
        Test Organism nodes
        :return:
        """
        df = Organism.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 600)
        self.assertEqual(properties, 16)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.name.nunique(), entries)

        self.assertNotEqual(df.ncbi_taxonomy.nunique(), entries)
        self.assertEqual(df.ncbi_taxonomy.nunique(), df.ncbi_taxonomy.count())

    def test_pathway(self):
        """
        Test Pathway nodes
        :return:
        """
        df = Pathway.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 250)
        self.assertEqual(properties, 7)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.name.nunique(), entries)

    def test_publication(self):
        """
        Test Publication nodes
        :return:
        """
        df = Publication.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 18500)
        self.assertEqual(properties, 10)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.pmid.nunique(), entries)
        self.assertNotEqual(df.year.nunique(), entries)

    def test_regulator(self):
        """
        Test Regulator nodes
        :return:
        """
        df = Regulator.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 16100)
        self.assertEqual(properties, 21)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.locus_tag.nunique(), entries)

        uniprot_accessions = df.uniprot_accession.nunique()
        self.assertNotEqual(uniprot_accessions, entries)
        self.assertEqual(uniprot_accessions, df.uniprot_accession.count())
        self.assertGreaterEqual(uniprot_accessions / entries, 0.65)

    def test_rfam(self):
        """
        Test RegulatoryFamily nodes
        :return:
        """
        df = RegulatoryFamily.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 210)
        self.assertEqual(properties, 9)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.name.nunique(), entries)
        self.assertNotEqual(df.mechanism.nunique(), entries)

    def test_interaction(self):
        """
        Test RegulatoryInteraction nodes
        :return:
        """
        df = RegulatoryInteraction.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 210000)
        self.assertEqual(properties, 12)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.interaction_hash.nunique(), entries)
        self.assertEqual(df.regulatory_effect.nunique(), 4)

        self.assertNotEqual(df.regulator.nunique(), entries)
        self.assertNotEqual(df.gene.nunique(), entries)

        self.assertLessEqual(df.effector.nunique(), entries)
        self.assertLessEqual(df.tfbs.nunique(), entries)
        self.assertLessEqual(df.organism.nunique(), entries)

    def test_source(self):
        """
        Test Source nodes
        :return:
        """
        df = Source.node_to_df()
        entries, properties = df.shape
        self.assertEqual(entries, 10)
        self.assertEqual(properties, 11)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.name.nunique(), entries)
        self.assertNotEqual(df.type.nunique(), entries)

    def test_tfbs(self):
        """
        Test TFBS nodes
        :return:
        """
        df = TFBS.node_to_df()
        entries, properties = df.shape
        self.assertGreaterEqual(entries, 83000)
        self.assertEqual(properties, 12)
        self.assertEqual(df.protrend_id.nunique(), entries)
        self.assertEqual(df.site_hash.nunique(), entries)

        self.assertNotEqual(df.organism.nunique(), entries)
        self.assertNotEqual(df.length.nunique(), entries)

        self.assertLessEqual(df.sequence.nunique(), entries)

        self.assertEqual(df.sequence.count(), entries)
        self.assertEqual(df.length.count(), entries)
        self.assertEqual(df.organism.count(), entries)

        self.assertGreaterEqual(df.start.count() / entries, 0.75)
        self.assertGreaterEqual(df.stop.count() / entries, 0.75)


if __name__ == '__main__':
    unittest.main()
