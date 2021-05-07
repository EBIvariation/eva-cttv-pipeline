import unittest

import json
import os
import requests
import xml.etree.ElementTree as ElementTree

from eva_cttv_pipeline.clinvar_xml_utils import ClinVarTrait
from eva_cttv_pipeline.evidence_string_generation import clinvar_to_evidence_strings
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT
from tests.eva_cttv_pipeline.evidence_string_generation import config


EFO_MAPPINGS = clinvar_to_evidence_strings.load_efo_mapping(config.efo_mapping_file)
GENE_MAPPINGS = CT.process_consequence_type_file(config.snp_2_gene_file)


class GetMappingsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.efo_mappings = EFO_MAPPINGS
        cls.gene_mappings = GENE_MAPPINGS

    def test_efo_mapping(self):
        self.assertEqual(len(self.efo_mappings), 9)

        self.assertEqual(self.efo_mappings['renal-hepatic-pancreatic dysplasia 2'][0],
                         ('http://www.orpha.net/ORDO/Orphanet_294415', 'Renal-hepatic-pancreatic dysplasia'))
        self.assertEqual(self.efo_mappings['frontotemporal dementia, ubiquitin-positive'][0],
                         ('http://www.orpha.net/ORDO/Orphanet_282', 'Frontotemporal dementia'))
        self.assertEqual(
            self.efo_mappings['3 beta-hydroxysteroid dehydrogenase deficiency'][0],
            ('http://www.orpha.net/ORDO/Orphanet_90791',
             'Congenital adrenal hyperplasia due to 3-beta-hydroxysteroid dehydrogenase deficiency'))

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
        orig_allele_origins = ['germline']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['germline']], converted_allele_origins)

    def test_just_somatic(self):
        orig_allele_origins = ['somatic']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['somatic']], converted_allele_origins)

    def test_mixed_germline(self):
        orig_allele_origins = ['germline', 'inherited', 'not applicable']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['germline', 'inherited', 'not applicable']], converted_allele_origins)

    def test_duplicate(self):
        orig_allele_origins = ['germline', 'germline']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['germline']], converted_allele_origins)
        orig_allele_origins = ['inherited', 'inherited', 'germline']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['germline', 'inherited']], converted_allele_origins)
        orig_allele_origins = ['somatic', 'somatic', 'somatic']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['somatic']], converted_allele_origins)

    def test_stringcase(self):
        orig_allele_origins = ['Germline']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['germline']], converted_allele_origins)
        orig_allele_origins = ['SOMATIC']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['somatic']], converted_allele_origins)

    def test_mixed(self):
        orig_allele_origins = ['germline', 'somatic']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['somatic'], ['germline']], converted_allele_origins)
        orig_allele_origins = ['somatic', 'inherited', 'not applicable']
        converted_allele_origins = clinvar_to_evidence_strings.convert_allele_origins(orig_allele_origins)
        self.assertListEqual([['somatic'], ['inherited', 'not applicable']], converted_allele_origins)


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
            []
        )


class GenerateEvidenceStringTest(unittest.TestCase):
    """Verifies that the evidence strings generated from a test ClinVar record matches the expectation."""

    def setUp(self):
        # This is required to display evidence string diffs (if any)
        self.maxDiff = None
        # Attributes for the evidence string generation
        self.clinvar_record = config.get_test_clinvar_record()
        self.disease_name = 'Rare congenital non-syndromic heart malformation'
        self.disease_source_id = 'C4017284'
        self.disease_mapped_efo_id = 'Orphanet_88991'
        self.consequence_attributes = GENE_MAPPINGS['14:67729209:A:G'][0]
        # Open Targets JSON schema
        schema_url = f'https://raw.githubusercontent.com/opentargets/json_schema/{config.OT_SCHEMA_VERSION}/' \
                     f'opentargets.json'
        self.ot_schema_contents = requests.get(schema_url).json()

    def test_genetics_evidence_string(self):
        """Verifies expected genetics evidence string generation."""
        evidence = clinvar_to_evidence_strings.generate_evidence_string(
            clinvar_record=self.clinvar_record,
            allele_origins=['germline'],
            disease_name=self.disease_name,
            disease_source_id=self.disease_source_id,
            disease_mapped_efo_id=self.disease_mapped_efo_id,
            consequence_attributes=self.consequence_attributes
        )
        # Check that the evidence string validates against schema
        clinvar_to_evidence_strings.validate_evidence_string(evidence, self.ot_schema_contents)
        # Check that the evidence string contents are as expected
        evidence_string = json.dumps(evidence, sort_keys=True, indent=2)
        expected_evidence_string = open(config.expected_genetics_evidence_string).read()
        self.assertEqual(evidence_string, expected_evidence_string)

    def test_somatic_evidence_string(self):
        """Verifies expected somatic evidence string generation."""
        evidence = clinvar_to_evidence_strings.generate_evidence_string(
            clinvar_record=self.clinvar_record,
            allele_origins=['somatic'],
            disease_name=self.disease_name,
            disease_source_id=self.disease_source_id,
            disease_mapped_efo_id=self.disease_mapped_efo_id,
            consequence_attributes=self.consequence_attributes
        )
        # Check that the evidence string validates against schema
        clinvar_to_evidence_strings.validate_evidence_string(evidence, self.ot_schema_contents)
        # Check that the evidence string contents are as expected
        evidence_string = json.dumps(evidence, sort_keys=True, indent=2)
        expected_evidence_string = open(config.expected_somatic_evidence_string).read()
        self.assertEqual(evidence_string, expected_evidence_string)


class GroupDiseasesByMappingTest(unittest.TestCase):
    """Verifies behaviour of group_diseases_by_efo_mapping"""

    def setUp(self):
        self.string_to_efo_mappings = {
            'disease a': [('EFO_1', 'Term 1')],
            'disease b': [('EFO_1', 'Term 1')],
            'disease c': [('EFO_1', 'Term 1')],
            'disease d': [('EFO_2', 'Term 2'), ('EFO_3', 'Term 3')]
        }

    @staticmethod
    def get_trait(name, medgen):
        """Return ClinVarTrait with given name and medgen id"""
        # template for trait xml that allows name and medgen id to be inserted
        trait_xml_template = f'''
        <Trait ID="123" Type="Disease">
            <Name>
                <ElementValue Type="Preferred">{name}</ElementValue>
            </Name>
            <XRef ID="{medgen}" DB="MedGen" />
        </Trait>
        '''
        return ClinVarTrait(trait_xml=ElementTree.fromstring(trait_xml_template), clinvar_record=None)

    def test_multiple_diseases_mapped_to_one_term(self):
        # Diseases mapped to same EFO term should be grouped together,
        # and the lexicographically first one returned.
        trait_a = self.get_trait('Disease A', 'MedGen_A')
        trait_b = self.get_trait('Disease B', 'MedGen_B')
        trait_c = self.get_trait('Disease C', 'MedGen_C')
        expected_result = [('Disease A', 'MedGen_A', 'EFO_1')]
        result = clinvar_to_evidence_strings.group_diseases_by_efo_mapping(
            clinvar_record_traits=[trait_b, trait_c, trait_a],
            string_to_efo_mappings=self.string_to_efo_mappings,
        )
        self.assertEqual(result, expected_result)

    def test_one_disease_mapped_to_multiple_terms(self):
        # One disease mapped to multiple terms should be split.
        trait_d = self.get_trait('Disease D', 'MedGen_D')
        expected_result = [('Disease D', 'MedGen_D', 'EFO_2'), ('Disease D', 'MedGen_D', 'EFO_3')]
        result = clinvar_to_evidence_strings.group_diseases_by_efo_mapping(
            clinvar_record_traits=[trait_d],
            string_to_efo_mappings=self.string_to_efo_mappings,
        )
        self.assertEqual(result, expected_result)

    def test_disease_without_mapping(self):
        # Diseases without mappings should be skipped.
        trait_e = self.get_trait('Disease E', 'MedGen_E')
        expected_result = []
        result = clinvar_to_evidence_strings.group_diseases_by_efo_mapping(
            clinvar_record_traits=[trait_e],
            string_to_efo_mappings=self.string_to_efo_mappings,
        )
        self.assertEqual(result, expected_result)
