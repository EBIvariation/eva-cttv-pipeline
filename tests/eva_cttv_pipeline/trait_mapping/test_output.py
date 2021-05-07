import csv
import tempfile
import unittest

import eva_cttv_pipeline.trait_mapping.output as output
from eva_cttv_pipeline.trait_mapping.oxo import OxOMapping, OxOResult
from eva_cttv_pipeline.trait_mapping.trait import OntologyEntry, Trait
import eva_cttv_pipeline.trait_mapping.zooma as zooma


class TestOutputTraitMapping(unittest.TestCase):
    def test_output_trait_mapping(self):
        tempfile_path = tempfile.mkstemp()[1]
        with open(tempfile_path, "w", newline='') as mapping_file:
            mapping_writer = csv.writer(mapping_file, delimiter="\t")
            mapping_writer.writerow(["#clinvar_trait_name", "uri", "label"])

            test_trait = Trait('aprt deficiency, japanese type', 11)

            # Normally a set, but changed to a list for predictable output order in test
            test_trait.finished_mapping_set = [
                OntologyEntry('http://www.orpha.net/ORDO/Orphanet_976',
                              'Adenine phosphoribosyltransferase deficiency'),
                OntologyEntry('http://www.orpha.net/ORDO/Orphanet_977',
                              'Adenine phosphoribosyltransferase deficiency type A')
            ]

            output.output_trait_mapping(test_trait, mapping_writer)

        with open(tempfile_path, "rt", newline='') as mapping_file:
            mapping_reader = csv.reader(mapping_file, delimiter="\t")
            next(mapping_reader)
            self.assertEqual(['aprt deficiency, japanese type',
                              'http://www.orpha.net/ORDO/Orphanet_976',
                              'Adenine phosphoribosyltransferase deficiency'],
                             next(mapping_reader))

            self.assertEqual(['aprt deficiency, japanese type',
                              'http://www.orpha.net/ORDO/Orphanet_977',
                              'Adenine phosphoribosyltransferase deficiency type A'],
                             next(mapping_reader))


class TestGetMappingsForCuration(unittest.TestCase):
    def test_get_non_efo_mapping(self):
        """If mapping is not in EFO, its `is_current` flag should *not* be checked, and the mapping
        *should* be selected for curation."""
        test_zooma_result = zooma.ZoomaResult(['http://purl.obolibrary.org/obo/HP_0001892'],
                                              'abnormal bleeding', 'HIGH', 'eva-clinvar')
        mapping = test_zooma_result.mapping_list[0]
        mapping.confidence = zooma.ZoomaConfidence.HIGH
        mapping.in_efo = False
        mapping.is_current = False
        mapping.ontology_label = ""
        mapping.source = 'eva-clinvar'
        mapping.uri = 'http://purl.obolibrary.org/obo/HP_0000483'
        self.assertEqual([mapping], output.get_mappings_for_curation([test_zooma_result]))

    def test_get_obsolete_efo_mapping(self):
        """If mapping is in EFO, but is not current, it *should not* be selected for curation."""
        test_zooma_result = zooma.ZoomaResult(['http://www.orpha.net/ORDO/Orphanet_976'],
                                              'Adenine phosphoribosyltransferase deficiency',
                                              'HIGH', 'eva-clinvar')
        mapping = test_zooma_result.mapping_list[0]
        mapping.confidence = zooma.ZoomaConfidence.HIGH
        mapping.in_efo = True
        mapping.is_current = False
        mapping.ontology_label = "Adenine phosphoribosyltransferase deficiency"
        mapping.source = 'eva-clinvar'
        mapping.uri = 'http://www.orpha.net/ORDO/Orphanet_976'
        self.assertEqual([], output.get_mappings_for_curation([test_zooma_result]))

    def test_get_current_efo_mapping(self):
        """If mapping is in EFO and is current, is *should* be selected for curation."""
        test_zooma_result = zooma.ZoomaResult(['http://purl.obolibrary.org/obo/MONDO_0008091'],
                                              'Abnormal neutrophil chemotactic response',
                                              'MEDIUM', 'eva-clinvar')
        mapping = test_zooma_result.mapping_list[0]
        mapping.confidence = zooma.ZoomaConfidence.HIGH
        mapping.in_efo = True
        mapping.is_current = True
        mapping.ontology_label = "Abnormal neutrophil chemotactic response"
        mapping.source = 'eva-clinvar'
        mapping.uri = 'http://purl.obolibrary.org/obo/MONDO_0008091'
        self.assertEqual([mapping], output.get_mappings_for_curation([test_zooma_result]))


class TestOutputForCuration(unittest.TestCase):
    def test_output_for_curation(self):
        tempfile_path = tempfile.mkstemp()[1]
        with open(tempfile_path, "wt") as curation_file:
            curation_writer = csv.writer(curation_file, delimiter="\t")

            test_trait = Trait("transitional cell carcinoma of the bladder", 276)

            test_oxo_result = OxOResult("HP:0006740", "Transitional cell carcinoma of the bladder",
                                        "HP:0006740")
            test_oxo_mapping = OxOMapping("bladder transitional cell carcinoma", "EFO:0006544", 2,
                                          "HP:0006740")
            test_oxo_mapping.in_efo = test_oxo_mapping.is_current = True
            test_oxo_mapping.ontology_label = "bladder transitional cell carcinoma"
            test_oxo_result.mapping_list = [test_oxo_mapping]

            test_trait.oxo_result_list = [test_oxo_result]

            output.output_for_curation(test_trait, curation_writer)

        with open(tempfile_path, "rt") as curation_file:
            curation_reader = csv.reader(curation_file, delimiter="\t")
            expected_record = [
                "transitional cell carcinoma of the bladder", "276", '',
                "http://www.ebi.ac.uk/efo/EFO_0006544|bladder transitional cell carcinoma|2|HP:0006740|EFO_CURRENT"
            ]
            self.assertEqual(expected_record, next(curation_reader))


if __name__ == '__main__':
    unittest.main()
