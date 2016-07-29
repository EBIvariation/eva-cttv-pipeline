from datetime import datetime
import json
import os
import unittest

from eva_cttv_pipeline import clinvar
from eva_cttv_pipeline import consequence_type as CT
from eva_cttv_pipeline import utilities

from tests import test_clinvar_to_evidence_strings
import tests.test_config as test_config


class TestClinvarRecord(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_clinvar_record = get_test_record()
        cls.rcv_to_rs, cls.rcv_to_nsv = \
            clinvar.get_rcv_to_rsnsv_mapping(test_config.variant_summary_file)
        cls.consequence_type_dict = CT.process_consequence_type_file(test_config.snp_2_gene_file)

    def test_gene_id(self):
        self.assertEqual(self.test_clinvar_record.gene_id, "NM_174878")

    def test_ensembl_id(self):
        self.assertEqual(self.test_clinvar_record.ensembl_id, "ENSG00000163646")

    def test_date(self):
        self.assertEqual(self.test_clinvar_record.date,
                         datetime.fromtimestamp(1435359600000/1000).isoformat())

    def test_score(self):
        self.assertEqual(self.test_clinvar_record.score, 1)

    def test_acc(self):
        self.assertEqual(self.test_clinvar_record.accession, "RCV000002127")

    def test_traits(self):
        self.assertEqual(self.test_clinvar_record.traits, [['Usher syndrome, type 3',
                                                            'Usher Syndrome, Type III',
                                                            'USHER SYNDROME, TYPE IIIA',
                                                            'Usher syndrome, type 3A']])

    def test_trait_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.trait_pubmed_refs, [[21697857]])

    def test_observed_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.observed_pubmed_refs, [11524702, 12145752])

    def test_measure_set_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.measure_set_pubmed_refs, [])

    def test_hgvs(self):
        self.assertEqual(self.test_clinvar_record.hgvs,
                         ['NM_174878.2:c.528T>G',
                          'NM_001256819.1:c.*142T>G',
                          'NM_052995.2:c.300T>G',
                          'NM_001195794.1:c.567T>G',
                          'NG_009168.1:g.49893T>G',
                          'NC_000003.12:g.150928107A>C',
                          'NC_000003.11:g.150645894A>C',
                          'NR_046380.2:n.1009T>G',
                          'NR_046380.1:n.1010T>G',
                          'p.Tyr176X',
                          'NP_443721.1:p.Tyr100Ter',
                          'NP_777367.1:p.Tyr176Ter',
                          'NP_001182723.1:p.Tyr189Ter'])

    def test_clinical_significance(self):
        self.assertEqual(self.test_clinvar_record.clinical_significance, "pathogenic")

    def test_get_rs(self):
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_rs(self.rcv_to_rs), "rs28940313")
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_rs({}), None)

    def test_get_nsv(self):
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_nsv(self.rcv_to_nsv), None)
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_nsv({"RCV000002127": "nsv123test"}),
                         "nsv123test")

    def test___get_main_consequence_types(self):
        test_consequence_type = CT.ConsequenceType(ensembl_gene_ids=["ENSG00000139988"],
                                                   so_names=["sequence_variant"])

        # print([so_name.__dict__ for so_name in self.consequence_type_dict["rs28940313"].so_terms])
        # print(self.rcv_to_rs)

        self.assertEqual(
            self.test_clinvar_record._ClinvarRecord__get_main_consequence_types(
                self.consequence_type_dict, self.rcv_to_rs),
            test_consequence_type)
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_main_consequence_types(
            {}, {}),
            None)

    def test_variant_type(self):
        self.assertEqual(self.test_clinvar_record.variant_type, "single nucleotide variant")

    def test_allele_origins(self):
        self.assertEqual(self.test_clinvar_record.allele_origins, ['germline'])


class TestGetRcvToRSNSVMapping(unittest.TestCase):
    variant_summary_file_path = os.path.join(os.path.dirname(__file__), 'resources',
                                             'variant_summary_2015-05_test_extract.txt')
    rcv_to_rs, rcv_to_nsv = clinvar.get_rcv_to_rsnsv_mapping(variant_summary_file_path)

    def test_rcv_to_rs(self):
        self.assertEqual(self.rcv_to_rs["RCV000000012"], "rs397704705")
        self.assertEqual(self.rcv_to_rs["RCV000000381"], "rs137854556")
        self.assertEqual(self.rcv_to_rs["RCV000000204"], "rs121965059")

    def test_rcv_to_nsv(self):
        self.assertEqual(self.rcv_to_nsv["RCV000004182"], "nsv1067860")
        self.assertEqual(self.rcv_to_nsv["RCV000004183"], "nsv1067861")
        self.assertEqual(self.rcv_to_nsv["RCV000004554"], "nsv1067916")


def get_test_record():
    test_clinvar_record_filepath = os.path.join(os.path.dirname(__file__), 'resources',
                                              'test_clinvar_record.json')
    with utilities.open_file(test_clinvar_record_filepath, "rt") as f:
        test_record_dict = json.load(f)
    test_record = clinvar.ClinvarRecord(test_clinvar_to_evidence_strings.MAPPINGS,
                                        test_record_dict)
    return test_record

