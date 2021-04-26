import unittest

import requests_mock

import eva_cttv_pipeline.trait_mapping.ols as ols
import tests.eva_cttv_pipeline.trait_mapping.resources.test_ols_data as test_ols_data


class TestGetTraitNames(unittest.TestCase):
    def test_get_ontology_label_from_ols(self):
        url = "http://www.orpha.net/ORDO/Orphanet_199318"
        ols_request_url = ols.build_ols_query(url)
        with requests_mock.mock() as m:
            m.get(ols_request_url, json=test_ols_data.TestGetTraitNamesData.orphanet_199318_ols_terms_json)
            self.assertEqual(ols.get_ontology_label_from_ols(url),
                             '15q13.3 microdeletion syndrome')


class TestIsCurrentAndInEfo(unittest.TestCase):
    def test_is_current_and_in_efo(self):

        with requests_mock.mock() as m:
            url = "https://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_425"
            m.get(url,
                  json=test_ols_data.TestIsCurrentAndInEfoData.orphanet_425_ols_efo_json)

            self.assertEqual(ols.is_current_and_in_efo("http://www.orpha.net/ORDO/Orphanet_425"),
                             True)


class TestIsInEfo(unittest.TestCase):
    def test_is_in_efo(self):
        with requests_mock.mock() as m:
            url = "https://www.ebi.ac.uk/ols/api/ontologies/efo/terms/http%253A%252F%252Fwww.orpha.net%252FORDO%252FOrphanet_425"
            m.get(url,
                  json=test_ols_data.TestIsInEfoData.orphanet_425_ols_efo_json)

            self.assertEqual(ols.is_in_efo("http://www.orpha.net/ORDO/Orphanet_425"),
                             True)
