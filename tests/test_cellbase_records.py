import os
import unittest

from eva_cttv_pipeline import cellbase_records, config


@unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                 "Skipping this test on Travis CI.")
class CellbaseRecordsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cb_records = cellbase_records.CellbaseRecords()

    def test_get_curr_result_list(self):
        curr_result_list = self.cb_records._CellbaseRecords__get_curr_result_list()
        self.assertEqual(len(curr_result_list), config.BATCH_SIZE)

        curr_response = self.cb_records._CellbaseRecords__get_curr_response()
        self.cb_records.skip = curr_response['numTotalResults']
        curr_result_list_test = self.cb_records._CellbaseRecords__get_curr_result_list()
        self.assertIsNone(curr_result_list_test)

        self.cb_records.skip = 999999
        curr_result_list = self.cb_records._CellbaseRecords__get_curr_result_list()
        self.assertIsNone(curr_result_list)

        self.cb_records.skip = 0
        curr_response = self.cb_records._CellbaseRecords__get_curr_response()
        len_to_expect = 20
        self.cb_records.skip = curr_response['numTotalResults'] - len_to_expect
        curr_result_list = self.cb_records._CellbaseRecords__get_curr_result_list()
        self.assertEqual(len(curr_result_list), len_to_expect)
