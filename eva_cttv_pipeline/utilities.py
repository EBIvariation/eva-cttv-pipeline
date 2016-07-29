import argparse
import errno
import gzip
import subprocess
import sys
import importlib
import importlib._bootstrap
import importlib.util
import os
import shutil

import eva_cttv_pipeline.config as config


def open_file(file_path, mode):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


def get_resource_file(package, resource):
    spec = importlib.util.find_spec(package)
    if spec is None:
        return None
    mod = (sys.modules.get(package) or importlib._bootstrap._load(spec))
    if mod is None or not hasattr(mod, '__file__'):
        return None

    parts = resource.split('/')
    parts.insert(0, os.path.dirname(mod.__file__))
    resource_name = os.path.join(*parts)
    return resource_name


def copy_and_overwrite(from_path, to_path):
    if os.path.exists(to_path):
        shutil.rmtree(to_path)
    shutil.copytree(from_path, to_path)


def copy_dir(src, dest):
    try:
        copy_and_overwrite(src, dest)
    except OSError as exception:
        # If the error was caused because the source wasn't a directory
        if exception.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % exception)


def change_json_refs(local_schema_dir):

    command = "find " + local_schema_dir + \
              " -type f -exec sed -i -e " \
              "\"s/https:\/\/raw.githubusercontent.com\/CTTV\/json_schema\/master/file:\/\/" + \
              local_schema_dir.replace("/", "\/") + "/g\" {} \;"
    print("Carrying out command:\n" + command)
    subprocess.check_output(command, shell=True)

    command = "find " + local_schema_dir + \
              " -type f -exec sed -i -e " \
              "\"s/https:\/\/raw.githubusercontent.com\/OpenTargets\/json_schema\/master/file:\/\/" + \
              local_schema_dir.replace("/", "\/") + "/g\" {} \;"
    print("Carrying out command:\n" + command)
    subprocess.check_output(command, shell=True)

    evidence_base_json = os.path.join(local_schema_dir, "src/evidence/base.json")
    evidence_base_json_temp = os.path.join(local_schema_dir, "src/evidence/base_temp.json")
    command = "grep -v '\"id\": \"base_evidence\"' " + evidence_base_json + " > " + \
              evidence_base_json_temp + \
              "; mv " + evidence_base_json_temp + " " + evidence_base_json
    print("Carrying out command:\n" + command)
    subprocess.check_output(command, shell=True)

    command = "grep -v '\"id\": \"#single_lit_reference\"' " + evidence_base_json + " > " + \
              evidence_base_json_temp + \
              "; mv " + evidence_base_json_temp + " " + evidence_base_json
    print("Carrying out command:\n" + command)
    subprocess.check_output(command, shell=True)

    command = "find " + local_schema_dir + " -type f -exec sed -i -e " + \
              "\"s/evidence\/base.json#base_evidence\/definitions\/single_lit_reference/evidence" + \
              "\/base.json#definitions\/single_lit_reference/g\" {} \;"
    # command = "sed -i -e \"s/evidence\/base.json#base_evidence\/definitions\/
    # single_lit_reference/evidence\/base.json#definitions\/
    # single_lit_reference/g\" " + evidence_base_json
    print("Carrying out command:\n" + command)
    subprocess.check_output(command, shell=True)

    command = "rm -rf " + local_schema_dir + ".git " + local_schema_dir + ".gitignore"
    # command = "sed -i -e \"s/evidence\/base.json#base_evidence\/definitions\/
    # single_lit_reference/evidence\/base.json#definitions\/
    # single_lit_reference/g\" " + evidence_base_json
    print("Carrying out command:\n" + command)
    subprocess.check_output(command, shell=True)


def create_local_schema():
    json_schema_dir = get_resource_file(__package__, "resources/json_schema")
    # local_schema_dir = str(os.path.join(str(Path(json_schema_dir).parent), "schema_copy"))
    local_schema_dir = str(os.path.join(os.path.dirname(os.path.dirname(json_schema_dir)),
                                        config.LOCAL_SCHEMA))

    ready_file = local_schema_dir + "/READY"
    if os.path.exists(ready_file):
        return

    os.makedirs(local_schema_dir, exist_ok=True)
    print("Copying json schema")
    print(json_schema_dir + " to " + local_schema_dir)
    copy_dir(json_schema_dir, local_schema_dir)

    change_json_refs(local_schema_dir)

    open(ready_file, 'a').close()


def check_for_local_schema():
    local_schema = get_resource_file(__package__, config.LOCAL_SCHEMA)
    if not os.path.exists(local_schema):
        create_local_schema()


class ArgParser:

    """
    Uses argparse module to parse command line arguments.
    Arguments are used in the pipeline, including input file paths, output path, path to files to
    specify EFO urls to either ignore or alter, and clinical significances that will be allowed to
    generate evidence strings.
    """

    def __init__(self, argv):
        usage = """
        *******************************************************************************************
        Task: generate CTTV evidence strings from ClinVar mongo
        *******************************************************************************************
        """
        parser = argparse.ArgumentParser(usage)

        parser.add_argument("--clinSig", dest="clinical_significance",
                            help="""Optional. String containing a comma-sparated list with the
                            clinical significances that will be allowed to generate
                            evidence-strings. By default all clinical significances will be
                            considered. Possible tags: 'unknown','untested','non-pathogenic',
                            'probable-non-pathogenic','probable-pathogenic','pathogenic',
                            'drug-response','drug response','histocompatibility','other','benign',
                            'protective','not provided','likely benign','confers sensitivity',
                            'uncertain significance','likely pathogenic',
                            'conflicting data from submitters','risk factor','association' """,
                            default="pathogenic,likely pathogenic")
        parser.add_argument("--ignore", dest="ignore_terms_file",
                            help="""Optional. String containing full path to a txt file containing
                            a list of term urls which will be ignored during batch processing """,
                            default=None)
        parser.add_argument("--adapt", dest="adapt_terms_file",
                            help="""Optional. String containing full path to a txt file containing
                            a list of invalid EFO urls which will be adapted to a general valid url
                             during batch processing """,
                            default=None)
        parser.add_argument("--out", dest="out",
                            help="""String containing the name of the file were
                            results will be stored.""",
                            required=True)

        parser.add_argument("-e", "--efoMapFile", dest="efo_mapping_file",
                            help="Path to file with trait name to url mappings", required=True)
        parser.add_argument("-g", "--snp2GeneFile", dest="snp_2_gene_file",
                            help="Path to file with RS id to ensembl gene ID and consequence "
                                 "mappings", required=True)
        parser.add_argument("-v", "--variantSummaryFile", dest="variant_summary_file",
                            help="Path to file with RS id to ensembl gene ID and consequence "
                                 "mappings", required=True)
        parser.add_argument("-j", dest="json_file", help="File containing Clinvar records json "
                                                         "strings in the format of documents in "
                                                         "Cellbase. One record per line.")

        args = parser.parse_args(args=argv[1:])

        self.clinical_significance = args.clinical_significance
        self.ignore_terms_file = args.ignore_terms_file
        self.adapt_terms_file = args.adapt_terms_file
        self.out = args.out
        self.efo_mapping_file = args.efo_mapping_file
        self.snp_2_gene_file = args.snp_2_gene_file
        self.variant_summary_file = args.variant_summary_file
        self.json_file = args.json_file


def check_dir_exists_create(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
