import unittest

from protrend.model import *
from protrend.utils import NeoDatabase
from protrend.utils.set_list import ProtrendSetList


class RelationalDatabaseTest(unittest.TestCase):
    """
    Testing the database status regarding random relationships
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

    def test_regulator_organisms(self):
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

    def test_regulator_genes(self):
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

    def test_regulator_tfbs(self):
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

    def test_regulator_effectors(self):
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
