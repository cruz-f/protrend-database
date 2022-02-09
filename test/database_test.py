import unittest

from protrend.model import *
from protrend.utils import NeoDatabase
from protrend.utils.set_list import ProtrendSetList


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

        self.assertLessEqual(df.regulator.nunique(), entries)
        self.assertLessEqual(df.gene.nunique(), entries)

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

    def test_interactions_organisms(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        ecoli_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(organism__exact=ecoli.protrend_id))
        ecoli_rels = ProtrendSetList(ecoli.regulatory_interaction.all())
        self.assertEqual(len(ecoli_interactions), len(ecoli_rels))

        backwards_rels = ProtrendSetList()
        for interaction in ecoli_interactions:
            backwards_rels.extend(interaction.data_organism.all())
        self.assertEqual(len(backwards_rels), 1)

    def test_interactions_regulators(self):
        """
        It tests the regulatory interactions and its relationships
        """
        fur = Regulator.nodes.get(uniprot_accession='P0A9A9')
        fur_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(regulator__exact=fur.protrend_id))
        fur_rels = ProtrendSetList(fur.regulatory_interaction.all())
        self.assertEqual(len(fur_interactions), len(fur_rels))

        backwards_rels = ProtrendSetList()
        for interaction in fur_interactions:
            backwards_rels.extend(interaction.data_regulator.all())
        self.assertEqual(len(backwards_rels), 1)

    def test_interactions_genes(self):
        """
        It tests the regulatory interactions and its relationships
        """
        cir_a = Gene.nodes.get(uniprot_accession='P17315')
        cir_a_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(gene__exact=cir_a.protrend_id))
        cir_a_rels = ProtrendSetList(cir_a.regulatory_interaction.all())
        self.assertEqual(len(cir_a_interactions), len(cir_a_rels))

        backwards_rels = ProtrendSetList()
        for interaction in cir_a_interactions:
            backwards_rels.extend(interaction.data_gene.all())
        self.assertEqual(len(backwards_rels), 1)

    def test_interactions_tfbs(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        ecoli_tfbs = TFBS.nodes.get(sequence='GATAATTGTTATCGTTTGC',
                                    organism=ecoli.protrend_id,
                                    stop=2246982)
        tfbs_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(tfbs__exact=ecoli_tfbs.protrend_id))
        tfbs_rels = ProtrendSetList(ecoli_tfbs.regulatory_interaction.all())
        self.assertEqual(len(tfbs_interactions), len(tfbs_rels))

        backwards_rels = ProtrendSetList()
        for interaction in tfbs_interactions:
            backwards_rels.extend(interaction.data_tfbs.all())
        self.assertEqual(len(backwards_rels), 1)

    def test_regulators_organisms(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        fur = Regulator.nodes.get(uniprot_accession='P0A9A9')
        fur_organisms = ProtrendSetList(fur.organism.all())
        self.assertEqual(len(fur_organisms), 1)
        self.assertEqual(fur_organisms[0].protrend_id, ecoli.protrend_id)

        fur_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(regulator__exact=fur.protrend_id))
        fur_interactions_organisms = set([interaction.organism for interaction in fur_interactions])
        self.assertEqual(len(fur_interactions_organisms), 1)

        backwards_rels = ProtrendSetList()
        for interaction in fur_interactions:
            backwards_rels.extend(interaction.data_organism.all())
        self.assertEqual(len(backwards_rels), 1)

    def test_regulators_genes(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        fur = Regulator.nodes.get(uniprot_accession='P0A9A9')
        fur_genes = ProtrendSetList(fur.gene.all())
        cir_a = Gene.nodes.get(uniprot_accession='P17315')
        cir_a_regulators = ProtrendSetList(cir_a.regulator.all())
        self.assertIn(cir_a.protrend_id, set([_gene.protrend_id for _gene in fur_genes]))
        self.assertIn(fur.protrend_id, set([_regulator.protrend_id for _regulator in cir_a_regulators]))

        fur_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(regulator__exact=fur.protrend_id,
                                                                              gene__exact=cir_a.protrend_id,
                                                                              organism__exact=ecoli.protrend_id))
        self.assertGreaterEqual(len(fur_interactions), 1)

    def test_regulators_tfbs(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        fur = Regulator.nodes.get(uniprot_accession='P0A9A9')
        fur_tfbs = ProtrendSetList(fur.tfbs.all())
        ecoli_tfbs = TFBS.nodes.get(sequence='GATAATTGTTATCGTTTGC',
                                    organism=ecoli.protrend_id,
                                    stop=2246982)
        ecoli_tfbs_regulators = ProtrendSetList(ecoli_tfbs.regulator.all())
        self.assertIn(ecoli_tfbs.protrend_id, set([_tfbs.protrend_id for _tfbs in fur_tfbs]))
        self.assertIn(fur.protrend_id, set([_regulator.protrend_id for _regulator in ecoli_tfbs_regulators]))

        fur_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(regulator__exact=fur.protrend_id,
                                                                              tfbs__exact=ecoli_tfbs.protrend_id,
                                                                              organism__exact=ecoli.protrend_id))
        self.assertGreaterEqual(len(fur_interactions), 1)

    def test_regulators_effectors(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        fur = Regulator.nodes.get(uniprot_accession='P0A9A9')
        fur_effectors = ProtrendSetList(fur.effector.all())
        iron = Effector.nodes.get(name='Iron ion, (Fe2+)')
        iron_regulators = ProtrendSetList(iron.regulator.all())
        self.assertIn(iron.protrend_id, set([_effector.protrend_id for _effector in fur_effectors]))
        self.assertIn(fur.protrend_id, set([_regulator.protrend_id for _regulator in iron_regulators]))

        fur_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(regulator__exact=fur.protrend_id,
                                                                              effector__exact=iron.protrend_id,
                                                                              organism__exact=ecoli.protrend_id))
        self.assertGreaterEqual(len(fur_interactions), 1)

    def test_genes_organisms(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        cir_a = Gene.nodes.get(uniprot_accession='P17315')
        cir_a_organisms = ProtrendSetList(cir_a.organism.all())
        self.assertEqual(len(cir_a_organisms), 1)
        self.assertEqual(cir_a_organisms[0].protrend_id, ecoli.protrend_id)

        cir_a_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(gene__exact=cir_a.protrend_id))
        cir_a_interactions_organisms = set([interaction.organism for interaction in cir_a_interactions])
        self.assertEqual(len(cir_a_interactions_organisms), 1)

        backwards_rels = ProtrendSetList()
        for interaction in cir_a_interactions:
            backwards_rels.extend(interaction.data_organism.all())
        self.assertEqual(len(backwards_rels), 1)

    def test_genes_tfbs(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        cir_a = Gene.nodes.get(uniprot_accession='P17315')
        cir_a_tfbs = ProtrendSetList(cir_a.tfbs.all())
        ecoli_tfbs = TFBS.nodes.get(sequence='GATAATTGTTATCGTTTGC',
                                    organism=ecoli.protrend_id,
                                    stop=2246982)
        ecoli_tfbs_genes = ProtrendSetList(ecoli_tfbs.gene.all())
        self.assertIn(ecoli_tfbs.protrend_id, set([_tfbs.protrend_id for _tfbs in cir_a_tfbs]))
        self.assertIn(cir_a.protrend_id, set([_gene.protrend_id for _gene in ecoli_tfbs_genes]))

        cir_a_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(gene__exact=cir_a.protrend_id,
                                                                                tfbs__exact=ecoli_tfbs.protrend_id,
                                                                                organism__exact=ecoli.protrend_id))
        self.assertGreaterEqual(len(cir_a_interactions), 1)

    def test_tfbs_organisms(self):
        """
        It tests the regulatory interactions and its relationships
        """
        ecoli = Organism.nodes.get(ncbi_taxonomy=511145)
        ecoli_tfbs = TFBS.nodes.get(sequence='GATAATTGTTATCGTTTGC',
                                    organism=ecoli.protrend_id,
                                    stop=2246982)
        ecoli_tfbs_organisms = ProtrendSetList(ecoli_tfbs.data_organism.all())
        self.assertEqual(len(ecoli_tfbs_organisms), 1)
        self.assertIn(ecoli.protrend_id, set([_organism.protrend_id for _organism in ecoli_tfbs_organisms]))

        ecoli_tfbs_interactions = ProtrendSetList(RegulatoryInteraction.nodes.filter(tfbs__exact=ecoli_tfbs.protrend_id,
                                                                                     organism__exact=ecoli.protrend_id))
        self.assertGreaterEqual(len(ecoli_tfbs_interactions), 1)


if __name__ == '__main__':
    unittest.main()
