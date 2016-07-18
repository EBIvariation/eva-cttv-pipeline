import os
import unittest
import shutil

from eva_cttv_pipeline import utilities


class GetResourceFileTest(unittest.TestCase):
    def test_get_resource_file_existent(self):
        self.assertTrue(utilities.get_resource_file("eva_cttv_pipeline", "resources/json_schema"))

    def test_get_resource_file_nonexistent(self):
        self.assertEqual(utilities.get_resource_file("not_a_real_package_39146", "not_a_real_file"),
                         None)


class CopyAndOverwriteTest(unittest.TestCase):
    def setUp(self):
        self.test_dir_a = os.path.join(os.path.dirname(__file__), "resources", "test_tmp_a")
        self.test_dir_b = os.path.join(os.path.dirname(__file__), "resources", "test_tmp_b")
        # self.test_dir_c = os.path.join(os.path.dirname(__file__), "resources", "test_tmp_c")
        os.makedirs(self.test_dir_a)
        os.makedirs(self.test_dir_b)
        self.test_file_a = os.path.join(self.test_dir_a, "test.txt")
        self.test_string = "this is a test string"
        with utilities.open_file(self.test_file_a, "wt") as f:
            f.write(self.test_string)

    def tearDown(self):
        shutil.rmtree(self.test_dir_a)
        shutil.rmtree(self.test_dir_b)
        # shutil.rmtree(self.test_dir_c)

    def test_existing(self):
        test_file_b = os.path.join(self.test_dir_b, "test.txt")
        with utilities.open_file(test_file_b, "wt") as f:
            f.write("hello world")
        utilities.copy_and_overwrite(self.test_dir_a, self.test_dir_b)
        with utilities.open_file(test_file_b, "rt") as f:
            contents = f.read()
        self.assertEqual(contents, self.test_string)


class CopyDirTest(unittest.TestCase):
    def setUp(self):
        self.test_dir_a = os.path.join(os.path.dirname(__file__), "resources", "test_tmp_a")
        self.test_dir_b = os.path.join(os.path.dirname(__file__), "resources", "test_tmp_b")
        self.test_dir_c = os.path.join(os.path.dirname(__file__), "resources", "test_tmp_c")
        os.makedirs(self.test_dir_a)
        os.makedirs(self.test_dir_b)
        self.test_file_a = os.path.join(self.test_dir_a, "test.txt")
        self.test_string = "this is a test string"
        with utilities.open_file(self.test_file_a, "wt") as f:
            f.write(self.test_string)

    def tearDown(self):
        shutil.rmtree(self.test_dir_a)
        shutil.rmtree(self.test_dir_b)

    def test_existing(self):
        test_file_b = os.path.join(self.test_dir_b, "test.txt")
        with utilities.open_file(test_file_b, "wt") as f:
            f.write("hello world")
        utilities.copy_dir(self.test_dir_a, self.test_dir_b)
        with utilities.open_file(test_file_b, "rt") as f:
            contents = f.read()
        self.assertEqual(contents, self.test_string)

    def test_nonexisting(self):
        utilities.copy_dir(self.test_dir_a, self.test_dir_c)
        with utilities.open_file(os.path.join(self.test_dir_c, "test.txt"), "rt") as f:
            contents = f.read()
        self.assertEqual(contents, self.test_string)
        shutil.rmtree(self.test_dir_c)


# TODO test change_json_refs, create_local_schema, check_for_local_schema


class ArgParserTest(unittest.TestCase):
    clin_sig = 'pathogenic,likely pathogenic'
    ignore = '/path/to/ignore/file'
    out = '/path/to/out/file'
    efo_map_file = '/path/to/efo/file'
    snp_2_gene_file = '/path/to/snp/to/gene/file'
    variant_summary_file = '/path/to/variant/summary/file'

    @classmethod
    def setUpClass(cls):
        argv = ['clinvar_to_evidence_strings.py', '--clinSig', cls.clin_sig, '--ignore',
                cls.ignore, '--out', cls.out, '-e', cls.efo_map_file, '-g', cls.snp_2_gene_file,
                '-v', cls.variant_summary_file]
        cls.argparser = utilities.ArgParser(argv)

    def test_clin_sig(self):
        self.assertEquals(self.argparser.clinical_significance, self.clin_sig)

    def test_ignore(self):
        self.assertEquals(self.argparser.ignore_terms_file, self.ignore)

    def test_out(self):
        self.assertEquals(self.argparser.out, self.out)

    def test_efo_map_file(self):
        self.assertEquals(self.argparser.efo_mapping_file, self.efo_map_file)

    def test_snp_2_gene_file(self):
        self.assertEquals(self.argparser.snp_2_gene_file, self.snp_2_gene_file)

    def test_variant_summary_file(self):
        self.assertEquals(self.argparser.variant_summary_file, self.variant_summary_file)


class CheckDirExistsCreateTest(unittest.TestCase):
    def test_create(self):
        directory = "./test_tmp"
        utilities.check_dir_exists_create(directory)
        self.assertTrue(os.path.exists(directory))
        os.rmdir(directory)
