import os
import unittest

from eva_cttv_pipeline import clinvar_to_evidence_strings
from tests import test_clinvar


def _get_mappings():
    efo_mapping_file = os.path.join(os.path.dirname(__file__), 'resources',
                                    'ClinVar_Traits_EFO_090915.xls')
    ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
    snp_2_gene_file = os.path.join(os.path.dirname(__file__), 'resources',
                                   'snp2gene_assignment_jul2016_extract.tsv.gz')
    variant_summary_file = os.path.join(os.path.dirname(__file__), 'resources',
                                        'variant_summary_2016-05_test_extract.txt')

    mappings = clinvar_to_evidence_strings.get_mappings(efo_mapping_file, ignore_file, None,
                                                        snp_2_gene_file, variant_summary_file)

    return mappings


MAPPINGS = _get_mappings()


class GetMappingsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mappings = MAPPINGS

    def test_efo_mapping(self):
        self.assertEqual(len(self.mappings.trait_2_efo), 3528)
        self.assertEqual(len(self.mappings.unavailable_efo_dict), 0)

        self.assertEqual(self.mappings.trait_2_efo["deafness, autosomal recessive 22"],
                         ['http://www.ebi.ac.uk/efo/EFO_0001063'])
        self.assertEqual(self.mappings.trait_2_efo["oculocutaneous albinism type 1b"],
                         ['http://www.orpha.net/ORDO/Orphanet_79434'])
        self.assertEqual(
            self.mappings.trait_2_efo["merosin deficient congenital muscular dystrophy"],
            ['http://www.orpha.net/ORDO/Orphanet_258'])

    def test_consequence_type_dict(self):
        self.assertEqual(len(self.mappings.consequence_type_dict), 54)

        self.assertTrue("rs724159824" in self.mappings.consequence_type_dict)
        self.assertTrue("rs34296458" in self.mappings.consequence_type_dict)
        self.assertTrue("rs199476100" in self.mappings.consequence_type_dict)
        self.assertTrue("rs80360485" in self.mappings.consequence_type_dict)

        self.assertFalse("rs0" in self.mappings.consequence_type_dict)
        self.assertFalse("rs5" in self.mappings.consequence_type_dict)
        self.assertFalse("rs9" in self.mappings.consequence_type_dict)

    def test_rcv_to_rs_nsv(self):
        self.assertEqual(len(self.mappings.rcv_to_rs), 18)
        self.assertEqual(len(self.mappings.rcv_to_nsv), 5)

        self.assertEqual(self.mappings.rcv_to_nsv["RCV000020147"], "nsv1067916")
        self.assertEqual(self.mappings.rcv_to_nsv["RCV000004182"], "nsv1067860")
        self.assertEqual(self.mappings.rcv_to_nsv["RCV000004183"], "nsv1067861")

        self.assertEqual(self.mappings.rcv_to_rs["RCV000061038"], "rs140870493")
        self.assertEqual(self.mappings.rcv_to_rs["RCV000038449"], "rs397517136")
        self.assertEqual(self.mappings.rcv_to_rs["RCV000126020"], "rs75686037")


class CreateTraitTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.trait = clinvar_to_evidence_strings.create_trait(9, ["Ciliary dyskinesia, primary, 7"],
                                                             MAPPINGS.trait_2_efo)

    def test_clinvar_trait_list(self):
        self.assertEqual(self.trait.clinvar_trait_list, ['Ciliary dyskinesia, primary, 7'])

    def test_efo_list(self):
        self.assertEqual(self.trait.efo_list, ['http://www.ebi.ac.uk/efo/EFO_0003900'])

    def test_return_none(self):
        none_trait = \
            clinvar_to_evidence_strings.create_trait(9, ["not a real trait"], MAPPINGS.trait_2_efo)
        self.assertIsNone(none_trait)


class SkipRecordTest(unittest.TestCase):

    def setUp(self):
        self.clinvar_record = test_clinvar.get_test_record()
        report = clinvar_to_evidence_strings.Report()
        # skip_record(clinvarRecord, cellbase_record, allowed_clinical_significance, counters)
        self.args = [self.clinvar_record, {"reference": "A", "alternate": "T"},
                     ["not provided"], report.counters]
        # allowed clin sig changed to just "non provided" to match that in the test record

    def test_return_false(self):
        self.assertFalse(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_ref_eq_alt(self):
        self.args[1] = {"reference": "A", "alternate": "A"}
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_rs_is_none(self):
        self.clinvar_record.rs = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))

    def test_con_type_is_none(self):
        self.clinvar_record.consequence_type = None
        self.assertTrue(clinvar_to_evidence_strings.skip_record(*self.args))


class LoadEfoMappingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ignore_file = os.path.join(os.path.dirname(__file__), 'resources', 'ignore_file.txt')
        efo_file = \
            os.path.join(os.path.dirname(__file__), 'resources', 'ClinVar_Traits_EFO_090915.xls')

        cls.trait_2_efo, cls.unavailable_efo = \
            clinvar_to_evidence_strings.load_efo_mapping(efo_file)
        cls.trait_2_efo_w_ignore, cls.unavailable_efo_w_ignore = \
            clinvar_to_evidence_strings.load_efo_mapping(efo_file, ignore_terms_file=ignore_file)

    def test_just_mapping_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo), 3819)

    def test_w_ignore_trait_2_efo(self):
        self.assertEqual(len(self.trait_2_efo_w_ignore), 3528)


class GetUnmappedUrlTest(unittest.TestCase):
    def test_orphanet(self):
        url = "http://www.orpha.net/ORDO/Orphanet_2670"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url),
                         "http://purl.bioontology.org/ORDO/Orphanet_2670")

    def test_hp(self):
        url = "http://purl.obolibrary.org/obo/HP_0000545"
        self.assertEqual(clinvar_to_evidence_strings.get_unmapped_url(url),
                         "http://purl.bioontology.org/obo/HP_0000545")

    def test_bad_url(self):
        #TODO currently the function exits execution with bad urls
        pass


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


class TestGetDefaultAllowedClincalSignificance(unittest.TestCase):
    def test_get_default_allowed_clincal_significance(self):
        correct_list = ['unknown', 'untested', 'non-pathogenic', 'probable-non-pathogenic',
         'probable-pathogenic', 'pathogenic', 'drug-response', 'drug response',
         'histocompatibility', 'other', 'benign', 'protective', 'not provided',
         'likely benign', 'confers sensitivity', 'uncertain significance',
         'likely pathogenic', 'conflicting data from submitters', 'risk factor',
         'association']
        self.assertEqual(clinvar_to_evidence_strings.get_default_allowed_clinical_significance(),
                         correct_list)
