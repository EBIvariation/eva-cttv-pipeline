import os
import unittest

from eva_cttv_pipeline import consequence_type


class ProcessGeneTest(unittest.TestCase):
    def test__process_gene(self):
        test_consequence_type_dict = {}
        test_rs_id = "rs121912888"
        test_ensembl_gene_id = "ENSG00000139219"
        test_so_name = "missense_variant"

        test_consequence_type = consequence_type.ConsequenceType(ensembl_gene_ids=[test_ensembl_gene_id],
                                                   so_names=[test_so_name])

        consequence_type.process_gene(test_consequence_type_dict, test_rs_id, test_ensembl_gene_id, test_so_name)

        self.assertEqual(test_consequence_type_dict["rs121912888"], test_consequence_type)


class ProcessConsequenceTypeFileTsvTest(unittest.TestCase):
    def test__process_consequence_type_file_tsv(self):
        test_consequence_type = consequence_type.ConsequenceType(ensembl_gene_ids=["ENSG00000021488"],
                                                   so_names=["missense_variant"])
        snp_2_gene_file_path = os.path.join(os.path.dirname(__file__), 'resources',
                                      'snp2gene_assignment_jul2016_extract.tsv')
        consequence_type_dict, one_rs_multiple_genes = \
            consequence_type.process_consequence_type_file_tsv(snp_2_gene_file_path)
        self.assertEqual(consequence_type_dict["rs121908485"], test_consequence_type)


class SoTermTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_so_term_a = consequence_type.SoTerm("stop_gained")
        cls.test_so_term_b = consequence_type.SoTerm("not_real_term")

    def test_accession(self):
        self.assertEqual(self.test_so_term_a.accession, "SO:0001587")
        self.assertIsNone(self.test_so_term_b.accession)

    def test_rank(self):
        self.assertEqual(self.test_so_term_a.rank, 3)
        self.assertEqual(self.test_so_term_b.rank, 34)


class ConsequenceTypeTest(unittest.TestCase):
    def setUp(self):
        self.test_ensembl_gene_ids_a = {"ENSG00000083093"}
        self.test_so_name_a = "intron_variant"
        self.test_so_name_b = "transcript_ablation"
        self.test_so_term_a = consequence_type.SoTerm(self.test_so_name_a)
        self.test_so_term_b = consequence_type.SoTerm(self.test_so_name_b)
        self.test_consequence_type_a = \
            consequence_type.ConsequenceType(ensembl_gene_ids=self.test_ensembl_gene_ids_a,
                               so_names={self.test_so_name_a})
        self.test_consequence_type_b = consequence_type.ConsequenceType()
        self.test_consequence_type_c = \
            consequence_type.ConsequenceType(ensembl_gene_ids=self.test_ensembl_gene_ids_a,
                               so_names={self.test_so_name_a, self.test_so_name_b})

    def test_ensembl_gene_ids(self):
        self.assertEqual(self.test_consequence_type_a.ensembl_gene_ids,
                         self.test_ensembl_gene_ids_a)
        new_id = "new_id"
        self.test_ensembl_gene_ids_a.add(new_id)
        self.test_consequence_type_a.ensembl_gene_ids.add(new_id)
        self.assertEqual(self.test_consequence_type_a.ensembl_gene_ids,
                         self.test_ensembl_gene_ids_a)

        new_ids = {"id1", "id2"}
        self.test_consequence_type_a.ensembl_gene_ids = new_ids
        self.assertEqual(self.test_consequence_type_a.ensembl_gene_ids, new_ids)

        self.assertEqual(self.test_consequence_type_b.ensembl_gene_ids, set())
        self.test_consequence_type_b.ensembl_gene_ids = new_ids
        self.assertEqual(self.test_consequence_type_b.ensembl_gene_ids, new_ids)

    def test_add_so_term(self):
        self.assertEqual(self.test_consequence_type_a.so_terms, {self.test_so_term_a})
        self.test_consequence_type_a.add_so_term(self.test_so_name_b)
        self.assertEqual(self.test_consequence_type_a.so_terms, {self.test_so_term_a,
                                                                 self.test_so_term_b})

        self.assertEqual(self.test_consequence_type_b.so_terms, set())
        self.test_consequence_type_b.add_so_term(self.test_so_name_b)
        self.assertEqual(self.test_consequence_type_b.so_terms, {self.test_so_term_b})

    def test_most_severe_so(self):
        self.assertEqual(self.test_consequence_type_a.most_severe_so, self.test_so_term_a)
        self.test_consequence_type_a.add_so_term(self.test_so_name_b)
        self.assertEqual(self.test_consequence_type_a.most_severe_so, self.test_so_term_b)

        self.assertEqual(self.test_consequence_type_c.most_severe_so, self.test_so_term_b)
