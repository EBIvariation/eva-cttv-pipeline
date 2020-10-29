import unittest

import os

from eva_cttv_pipeline.evidence_string_generation import clinvar_to_evidence_strings
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT
from tests.evidence_string_generation import config


EFO_MAPPINGS = clinvar_to_evidence_strings.load_efo_mapping(config.efo_mapping_file)
GENE_MAPPINGS = CT.process_consequence_type_file(config.snp_2_gene_file)


class GetMappingsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.efo_mappings = EFO_MAPPINGS
        cls.gene_mappings = GENE_MAPPINGS

    def test_efo_mapping(self):
        self.assertEqual(len(self.efo_mappings), 5283)

        self.assertEqual(self.efo_mappings['renal-hepatic-pancreatic dysplasia 2'][0],
                         ('http://www.orpha.net/ORDO/Orphanet_294415', None))
        self.assertEqual(self.efo_mappings['frontotemporal dementia'][0],
                         ('http://purl.obolibrary.org/obo/HP_0000733', None))
        self.assertEqual(
            self.efo_mappings['3 beta-hydroxysteroid dehydrogenase deficiency'][0],
            ('http://www.orpha.net/ORDO/Orphanet_90791', None))

        self.assertEqual(
            self.efo_mappings['coronary artery disease/myocardial infarction'],
            [('http://www.ebi.ac.uk/efo/EFO_0000612', 'myocardial infarction'),
             ('http://www.ebi.ac.uk/efo/EFO_0001645', 'coronary heart disease')])

    def test_consequence_type_dict(self):
        self.assertEqual(len(self.gene_mappings), 21)

        self.assertTrue('14:67727191:G:A' in self.gene_mappings)
        self.assertTrue('14:67727197:C:T' in self.gene_mappings)
        self.assertTrue('14:67729179:T:C' in self.gene_mappings)
        self.assertTrue('14:67729307:CG:C' in self.gene_mappings)

        self.assertFalse('rs0' in self.gene_mappings)
        self.assertFalse('rs5' in self.gene_mappings)
        self.assertFalse('rs9' in self.gene_mappings)


class GetTermsFromFileTest(unittest.TestCase):
    #TODO do the same for adapt terms file?
    @classmethod
    def setUpClass(cls):
        ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
        cls.ignore_terms = clinvar_to_evidence_strings.get_terms_from_file(ignore_file)

    def test_with_file(self):
        self.assertEqual(len(self.ignore_terms), 218)
        self.assertEqual(self.ignore_terms[0], "http://purl.obolibrary.org/obo/HP_0011677")
        self.assertEqual(self.ignore_terms[-1], "http://www.orpha.net/ORDO/Orphanet_120795")

    def test_no_file(self):
        self.assertEqual(clinvar_to_evidence_strings.get_terms_from_file(None), [])


class TestConvertAlleleOrigins(unittest.TestCase):
    def test_just_germline(self):
        orig_allele_origins = ["germline"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual(["germline"], converted_allele_origins)

    def test_just_somatic(self):
        orig_allele_origins = ["somatic"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual(["somatic"], converted_allele_origins)

    def test_just_tested_inconclusive(self):
        orig_allele_origins = ["tested-inconclusive"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([], converted_allele_origins)

    def test_just_other_germline(self):
        orig_allele_origins_list = [["unknown"],
                                    ["inherited"],
                                    ["maternal"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
            self.assertListEqual(["germline"], converted_allele_origins)

    def test_nonsense(self):
        orig_allele_origins_list = [["fgdsgfgs"],
                                    ["notarealorigin"],
                                    ["134312432:dasdfd"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
            self.assertListEqual([], converted_allele_origins)
        orig_allele_origins = ["fgdsgfgs", "germline"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual(["germline"], converted_allele_origins)

    def test_mixed_germline(self):
        orig_allele_origins_list = [["germline", "de novo"],
                                    ["germline", "inherited", "not applicable"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
            self.assertListEqual(["germline"], converted_allele_origins)

    def test_duplicate(self):
        orig_allele_origins = ["germline", "germline"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
            orig_allele_origins)
        self.assertListEqual(["germline"], converted_allele_origins)
        orig_allele_origins = ["inherited", "inherited", "germline"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
            orig_allele_origins)
        self.assertListEqual(["germline"], converted_allele_origins)
        orig_allele_origins = ["somatic", "somatic", "somatic"]
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
            orig_allele_origins)
        self.assertListEqual(["somatic"], converted_allele_origins)


    def test_stringcase(self):
        orig_allele_origins_list = [["Germline"],
                               ["InHerIted"],
                               ["UNKNOWN"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
            self.assertListEqual(["germline"], converted_allele_origins)
        orig_allele_origins_list = [["Somatic"],
                                    ["SOMATIC"],
                                    ["sOMatIc"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
                orig_allele_origins)
            self.assertListEqual(["somatic"], converted_allele_origins)

    def test_mixed(self):
        orig_allele_origins_list = [["germline", "somatic"],
                                    ["somatic", "inherited", "not applicable"],
                                    ["somatic", "unknown"]]
        for orig_allele_origins in orig_allele_origins_list:
            converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(
                orig_allele_origins)
            self.assertListEqual(["somatic", "germline"], converted_allele_origins)


class TestGetConsequenceTypes(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # A single example ClinVar record
        cls.test_crm = config.get_test_clinvar_record().measure
        # Example result from the gene & functional consequence mapping pipeline for several variants
        cls.consequence_type_dict = GENE_MAPPINGS

    def test_get_consequence_types(self):
        self.assertEqual(
            clinvar_to_evidence_strings.get_consequence_types(self.test_crm, self.consequence_type_dict)[0],
            CT.ConsequenceType('ENSG00000139988', CT.SoTerm('missense_variant')),
            ''
        )
        self.assertEqual(
            clinvar_to_evidence_strings.get_consequence_types(self.test_crm, {}),
            [None]
        )
