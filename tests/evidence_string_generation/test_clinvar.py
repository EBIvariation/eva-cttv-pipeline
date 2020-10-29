from datetime import datetime
import unittest


from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT
from tests.evidence_string_generation import config


class TestClinvarRecord(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_clinvar_record = config.get_test_clinvar_record()

    def test_date(self):
        """Check that the last updated date of the referenceClinVarAssertion is loaded correctly"""
        self.assertEqual(self.test_clinvar_record.date, '2020-09-16')

    def test_score(self):
        self.assertEqual(self.test_clinvar_record.score, 2)

    def test_review_status(self):
        self.assertEqual(self.test_clinvar_record.review_status, 'criteria provided, multiple submitters, no conflicts')

    def test_acc(self):
        self.assertEqual(self.test_clinvar_record.accession, 'RCV000002127')

    def test_traits(self):
        self.assertEqual(self.test_clinvar_record.traits[0].name, 'Leber congenital amaurosis 13')

    def test_trait_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.traits[0].pubmed_refs, [20301475, 20301590, 30285347])

    def test_observed_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.observed_pubmed_refs, [15258582, 15322982])

    def test_clinical_significance(self):
        self.assertEqual(self.test_clinvar_record.clinical_significance_list, ['likely pathogenic', 'pathogenic'])

    def test_allele_origins(self):
        self.assertEqual(self.test_clinvar_record.allele_origins, {'germline', 'inherited', 'unknown'})


class TestClinvarRecordMeasure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_crm = config.get_test_clinvar_record().measure
        cls.consequence_type_dict = CT.process_consequence_type_file(config.snp_2_gene_file)

    def test_hgvs(self):
        self.assertEqual(self.test_crm.hgvs,
                         ['NM_152443.3:c.677A>G',
                          'NG_008321.1:g.32324A>G',
                          'NC_000014.9:g.67729209A>G',
                          'NC_000014.8:g.68195926A>G',
                          'NM_152443.2:c.677A>G',
                          'Q96NR8:p.Tyr226Cys',
                          'NP_689656.2:p.Tyr226Cys'])

    def test_rs(self):
        self.assertEqual(self.test_crm.rs_id, 'rs28940313')

    def test_nsv(self):
        self.assertEqual(self.test_crm.nsv_id, None)

    def test_variant_type(self):
        self.assertEqual(self.test_crm.variant_type, 'single nucleotide variant')

    def test_measure_set_pubmed_refs(self):
        self.assertEqual(self.test_crm.pubmed_refs, [])
