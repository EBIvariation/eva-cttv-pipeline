import unittest

from eva_cttv_pipeline import file_utils


class GetResourceFileTest(unittest.TestCase):
    def test_get_resource_file_existent(self):
        self.assertTrue(file_utils.get_resource_file("eva_cttv_pipeline", "resources/json_schema"))

    def test_get_resource_file_nonexistent(self):
        self.assertEqual(file_utils.get_resource_file("not_a_real_package_39146", "not_a_real_file"),
                         None)
