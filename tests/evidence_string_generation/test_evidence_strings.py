from datetime import datetime
import gzip
import json
import os.path
from types import SimpleNamespace
import unittest

from eva_cttv_pipeline.evidence_string_generation import clinvar
from eva_cttv_pipeline.evidence_string_generation import clinvar_to_evidence_strings
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT
from eva_cttv_pipeline.evidence_string_generation import evidence_strings
from tests.evidence_string_generation import test_clinvar_to_evidence_strings
from tests.evidence_string_generation import config


def get_input_data_for_evidence_string_generation():
    """Prepares mock input data necessary for the evidence string generation."""
    clinvar_record = clinvar.ClinvarRecord(json.load(open(config.test_clinvar_record_file)))
    report = clinvar_to_evidence_strings.Report()

    trait = SimpleNamespace()
    trait.trait_counter = 0
    trait.clinvar_name = 'ClinVar name of a mock trait'
    trait.ontology_id = 'http://www.orpha.net/ORDO/Orphanet_88991'
    trait.ontology_label = None

    consequence_type = test_clinvar_to_evidence_strings.MAPPINGS.consequence_type_dict['14:67729209:A:G'][0]
    return clinvar_record, clinvar_record.measures[0], report, trait, consequence_type


class GenerateEvidenceStringTest(unittest.TestCase):
    """Verifies that the evidence strings generated from a test ClinVar record matches the expectation."""

    def setUp(self):
        self.maxDiff = None  # this is required to display evidence string diffs (if any)
        self.test_args = get_input_data_for_evidence_string_generation()

    def test_genetics_evidence_string(self):
        """Verifies expected genetics evidence string generation."""
        evidence_string = json.dumps(evidence_strings.CTTVGeneticsEvidenceString(*self.test_args),
                                     sort_keys=True, indent=2)
        expected_evidence_string = open(config.expected_genetics_evidence_string).read()
        self.assertEqual(evidence_string, expected_evidence_string)

    def test_somatic_evidence_string(self):
        """Verifies expected somatic evidence string generation."""
        evidence_string = json.dumps(evidence_strings.CTTVSomaticEvidenceString(*self.test_args),
                                     sort_keys=True, indent=2)
        expected_evidence_string = open(config.expected_somatic_evidence_string).read()
        self.assertEqual(evidence_string, expected_evidence_string)


class GeneticsEvidenceStringTest(unittest.TestCase):
    """Test the internal mechanics of the class which handles genetics evidence strings."""

    def setUp(self):
        self.test_args = get_input_data_for_evidence_string_generation()
        self.test_ges = evidence_strings.CTTVGeneticsEvidenceString(*self.test_args)
        self.ot_schema_contents = json.loads(gzip.open(config.open_targets_schema_gz).read().decode('utf-8'))

    def test_unique_association_field(self):
        uaf_1 = ("gene", "test_gene")
        uaf_2 = ("clinvarAccession", "test_clinvar")
        uaf_3 = ("alleleOrigin", "germline")
        uaf_4 = ("phenotype", "test_phenotype")
        uaf_5 = ("variant_id", "test_rs")

        self.test_ges.add_unique_association_field(*uaf_1)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_1[0]], uaf_1[1])
        self.test_ges.add_unique_association_field(*uaf_2)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_2[0]], uaf_2[1])

        self.test_ges.add_unique_association_field(*uaf_3)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_3[0]], uaf_3[1])
        self.test_ges.add_unique_association_field(*uaf_4)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_4[0]], uaf_4[1])
        self.test_ges.add_unique_association_field(*uaf_5)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_5[0]], uaf_5[1])

    def test_set_target(self):
        target = ("http://identifiers.org/ensembl/ENSG00000135486",
                  "http://identifiers.org/cttv.activity/predicted_damaging")
        self.test_ges._clear_target()
        self.test_ges.set_target(*target)
        self.assertEqual(self.test_ges['target']['id'], target[0])
        self.assertEqual(self.test_ges['target']['activity'], target[1])

    def test_disease(self):
        disease_id = "Ciliary dyskinesia, primary, 26"
        self.test_ges.disease_id = disease_id

    def test_evidence_codes(self):
        evidence_codes = ["http://purl.obolibrary.org/obo/ECO_0000205"]
        self.test_ges.evidence_codes = evidence_codes
        self.assertEqual(self.test_ges['evidence']['evidence_codes'], evidence_codes)
        self.assertEqual(self.test_ges.evidence_codes, evidence_codes)

    def test_top_level_literature(self):
        literature = ["http://europepmc.org/abstract/MED/20301537"]
        self.test_ges.top_level_literature = literature
        self.assertEqual(self.test_ges['literature']['references'],
                         [{"lit_id": literature_id} for literature_id in literature])
        self.assertEqual(self.test_ges.top_level_literature,
                         [{"lit_id": literature_id} for literature_id in literature])

    def test_db_xref_url(self):
        url = "http://identifiers.org/clinvar.record/RCV000128628"
        self.test_ges.db_xref_url = url
        self.assertEqual(
            self.test_ges['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'],
            url)
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'],
            url)
        self.assertEqual(self.test_ges.db_xref_url, url)

    def test_url(self):
        url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
        self.test_ges.url = url
        self.assertEqual(self.test_ges['evidence']['gene2variant']['urls'][0]['url'], url)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['urls'][0]['url'], url)
        self.assertEqual(self.test_ges.url, url)

    def test_gene_2_var_ev_codes(self):
        ev_codes = ['http://identifiers.org/eco/cttv_mapping_pipeline']
        self.test_ges.gene_2_var_ev_codes = ev_codes
        self.assertEqual(self.test_ges['evidence']['gene2variant']['evidence_codes'], ev_codes)
        self.assertEqual(self.test_ges.gene_2_var_ev_codes, ev_codes)

    def test_gene_2_var_func_consequence(self):
        functional_consequence = 'http://purl.obolibrary.org/obo/SO_0001583'
        self.test_ges.gene_2_var_func_consequence = functional_consequence
        self.assertEqual(self.test_ges['evidence']['gene2variant']['functional_consequence'],
                         functional_consequence)
        self.assertEqual(self.test_ges.gene_2_var_func_consequence, functional_consequence)

    def test_set_var_2_disease_literature_a(self):
        self.test_ges['evidence']['variant2disease']['provenance_type']['literature'] = {}

        literature_1 = "PMCID12345"
        self.test_ges.set_var_2_disease_literature([literature_1])
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'],
            [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ges.set_var_2_disease_literature(literature_list)
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'],
            [{"lit_id": literature_id} for literature_id in literature_list])

    def test_set_var_2_disease_literature_b(self):
        literature_1 = "PMCID12345"
        self.test_ges.set_var_2_disease_literature([literature_1])
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'],
            [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ges.set_var_2_disease_literature(literature_list)
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'],
            [{"lit_id": literature_id} for literature_id in literature_list])

    def test_association(self):
        self.test_ges.association = True
        self.assertTrue(self.test_ges['evidence']['gene2variant']['is_associated'])
        self.assertTrue(self.test_ges['evidence']['variant2disease']['is_associated'])
        self.assertTrue(self.test_ges.association)

        self.test_ges.association = False
        self.assertFalse(self.test_ges['evidence']['gene2variant']['is_associated'])
        self.assertFalse(self.test_ges['evidence']['variant2disease']['is_associated'])
        self.assertFalse(self.test_ges.association)

    def test_set_variant(self):
        test_id = "http://identifiers.org/dbsnp/rs193922494"
        test_type = "snp single"
        self.test_ges._clear_variant()
        self.test_ges.set_variant(test_id, test_type)
        self.assertEqual(self.test_ges['variant']['id'], test_id)
        self.assertEqual(self.test_ges['variant']['type'], test_type)

    def test_unique_reference(self):
        unique_reference = "http://europepmc.org/abstract/MED/0"
        self.test_ges.unique_reference = unique_reference
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['unique_experiment_reference'],
            unique_reference)
        self.assertEqual(self.test_ges.unique_reference, unique_reference)

    def test_date(self):
        date_string = datetime.fromtimestamp(1412982000000 / 1000).isoformat()
        self.test_ges.date = date_string
        self.assertEqual(self.test_ges['evidence']['gene2variant']['date_asserted'], date_string)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['date_asserted'],
                         date_string)
        self.assertEqual(self.test_ges.date, date_string)

    def test_validate(self):
        test_args = get_input_data_for_evidence_string_generation()
        test_evidence_string = evidence_strings.CTTVGeneticsEvidenceString(*test_args)
        self.assertTrue(test_evidence_string.validate(self.ot_schema_contents))


class SomaticEvidenceStringTest(unittest.TestCase):
    """Test the internal mechanics of the class which handles somatic evidence strings."""

    @classmethod
    def setUpClass(cls):
        cls.consequence_type_dict = CT.process_consequence_type_file(config.snp_2_gene_file)

    def setUp(self):
        test_args = get_input_data_for_evidence_string_generation()
        self.test_ses = evidence_strings.CTTVSomaticEvidenceString(*test_args)
        ot_schema_path = os.path.join(
            os.path.dirname(__file__), 'resources', 'opentargets.1.6.6.json.gz')
        self.ot_schema_contents = json.loads(gzip.open(ot_schema_path).read().decode('utf-8'))

    def test_db_xref_url(self):
        url = "http://identifiers.org/clinvar.record/RCV000128628"
        self.test_ses.db_xref_url = url
        self.assertEqual(self.test_ses['evidence']['provenance_type']['database']['dbxref']['url'],
                         url)
        self.assertEqual(self.test_ses.db_xref_url, url)

    def test_url(self):
        url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
        self.test_ses.url = url
        self.assertEqual(self.test_ses['evidence']['urls'][0]['url'], url)
        self.assertEqual(self.test_ses.url, url)

    def test_evidence_literature(self):
        literature_1 = "PMCID12345"
        self.test_ses.evidence_literature = [literature_1]
        self.assertEqual(self.test_ses['evidence']['provenance_type']['literature']['references'],
                         [{"lit_id": literature_1}])
        self.assertEqual(self.test_ses.evidence_literature, [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ses.evidence_literature = literature_list
        self.assertEqual(self.test_ses['evidence']['provenance_type']['literature']['references'],
                         [{"lit_id": literature_id} for literature_id in literature_list])
        self.assertEqual(self.test_ses.evidence_literature,
                         [{"lit_id": literature_id} for literature_id in literature_list])

    def test_association(self):
        self.test_ses.association = True
        self.assertTrue(self.test_ses['evidence']['is_associated'])
        self.assertTrue(self.test_ses.association)

        self.test_ses.association = False
        self.assertFalse(self.test_ses['evidence']['is_associated'])
        self.assertFalse(self.test_ses.association)

    def test_date(self):
        date_string = datetime.fromtimestamp(1412982000000 / 1000).isoformat()
        self.test_ses.date = date_string
        self.assertEqual(self.test_ses['evidence']['date_asserted'], date_string)
        self.assertEqual(self.test_ses.date, date_string)

    def test_add_known_mutations(self):
        functional_consequence = "http://purl.obolibrary.org/obo/SO_0001791"
        preferred_name = "exon_variant"
        self.test_ses._clear_known_mutations()
        self.test_ses.add_known_mutation(functional_consequence, preferred_name)
        self.assertEqual(
            self.test_ses['evidence']['known_mutations'],
            [{'functional_consequence': functional_consequence, 'preferred_name': preferred_name}])

    def test_set_known_mutations(self):
        test_consequence_type = CT.ConsequenceType("ENSG00000008710",
                                                   CT.SoTerm("3_prime_UTR_variant"))
        self.test_ses._clear_known_mutations()
        self.test_ses.set_known_mutations(test_consequence_type.so_term)
        self.assertEqual(
            self.test_ses['evidence']['known_mutations'],
            [{'functional_consequence': 'http://purl.obolibrary.org/obo/SO_0001624',
              'preferred_name': '3_prime_UTR_variant'}])

    def test_validate(self):
        test_args = get_input_data_for_evidence_string_generation()
        test_evidence_string = evidence_strings.CTTVSomaticEvidenceString(*test_args)
        self.assertTrue(test_evidence_string.validate(self.ot_schema_contents))


class GetOpenTargetsVariantTypeTest(unittest.TestCase):
    """Test determination of a variant type (as defined in the Open Targets schema) based on allele sequences."""
    @classmethod
    def setUpClass(cls):

        crm = SimpleNamespace()
        crm.vcf_ref = "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG"
        crm.vcf_alt = "C"
        record_single_a = (crm, "snp single")
        crm = SimpleNamespace()
        crm.vcf_ref = "A"
        crm.vcf_alt = "C"
        record_single_b = (crm, "snp single")
        record_single_c = SimpleNamespace()
        record_single_c.vcf_ref = "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG"
        record_single_c.vcf_alt = "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG"
        record_single_c = (crm, "snp single")

        cls.test_records_singles = [record_single_a, record_single_b, record_single_c]

        crm = SimpleNamespace()
        crm.vcf_ref = "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        crm.vcf_alt = "C"
        record_structurals_a = (crm, "structural variant")
        crm = SimpleNamespace()
        crm.vcf_ref = "A"
        crm.vcf_alt = "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        record_structurals_b = (crm, "structural variant")
        record_single_c = SimpleNamespace()
        record_single_c.vcf_ref = "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        record_single_c.vcf_alt = "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        record_structurals_c = (crm, "structural variant")

        cls.test_records_structurals = \
            [record_structurals_a, record_structurals_b, record_structurals_c]

    def test_get_cttv_variant_type_singles(self):
        for record in self.test_records_singles:
            self.assertEqual(evidence_strings.get_cttv_variant_type(record[0]), record[1])

    def test_get_cttv_variant_type_structurals(self):
        for record in self.test_records_structurals:
            self.assertEqual(evidence_strings.get_cttv_variant_type(record[0]), record[1])
