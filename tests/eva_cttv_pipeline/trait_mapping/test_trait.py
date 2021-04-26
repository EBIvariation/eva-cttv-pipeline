import unittest

import eva_cttv_pipeline.trait_mapping.trait as trait
import eva_cttv_pipeline.trait_mapping.zooma as zooma
import eva_cttv_pipeline.trait_mapping.oxo as oxo


class TestTrait(unittest.TestCase):
    def test_is_finished_true(self):
        test_trait = trait.Trait('aprt deficiency, japanese type', 1)
        test_trait.finished_mapping_set.add(trait.OntologyEntry('http://www.orpha.net/ORDO/Orphanet_976',
                                                                'Adenine phosphoribosyltransferase deficiency'))
        self.assertTrue(test_trait.is_finished)

    def test_is_finished_false(self):
        test_trait = trait.Trait('aprt deficiency, japanese type', 1)
        self.assertFalse(test_trait.is_finished)

    def test_process_zooma_result(self):
        test_trait = trait.Trait('aprt deficiency, japanese type', 1)

        test_zooma_result = zooma.ZoomaResult(['http://www.orpha.net/ORDO/Orphanet_976'],
                                              'Adenine phosphoribosyltransferase deficiency',
                                              'HIGH', 'eva-clinvar')
        entry = test_zooma_result.mapping_list[0]
        entry.confidence = zooma.ZoomaConfidence.HIGH
        entry.in_efo = True
        entry.is_current = True
        entry.ontology_label = "Adenine phosphoribosyltransferase deficiency"
        entry.source = 'eva-clinvar'
        entry.uri = 'http://www.orpha.net/ORDO/Orphanet_976'

        test_trait.zooma_result_list.append(test_zooma_result)

        test_trait.process_zooma_results()
        self.assertEqual(1, len(test_trait.finished_mapping_set))

    def test_process_oxo_mappings(self):
        test_trait = trait.Trait('congenital cystic disease of liver', 11)

        test_oxo_result = oxo.OxOResult('HP:0006706', 'Cystic liver disease', 'HP:0006706')

        test_oxo_mapping_1 = oxo.OxOMapping('Isolated polycystic liver disease', 'Orphanet:2924', 2,
                                            'HP:0006706')
        test_oxo_mapping_1.in_efo = True
        test_oxo_mapping_1.is_current = True

        test_oxo_mapping_2 = oxo.OxOMapping('cystic liver disease', 'EFO:1001505', 1, 'HP:0006706')
        test_oxo_mapping_2.in_efo = True
        test_oxo_mapping_2.is_current = True

        test_oxo_result.mapping_list = [test_oxo_mapping_1, test_oxo_mapping_2]

        test_trait.oxo_result_list.append(test_oxo_result)

        test_trait.process_oxo_mappings()

        self.assertEqual(1, len(test_trait.finished_mapping_set))




