import os

from eva_cttv_pipeline import clinvar_xml_utils

OT_SCHEMA_VERSION = "2.0.6"

test_dir = os.path.dirname(__file__)
efo_mapping_file = os.path.join(test_dir, 'resources', 'string_to_ontology_mappings.tsv')
snp_2_gene_file = os.path.join(test_dir, 'resources/snp2gene_extract.tsv')
test_clinvar_record_file = os.path.join(test_dir, 'resources/test_clinvar_record.xml.gz')
expected_genetics_evidence_string = os.path.join(test_dir, 'resources/expected_genetics_evidence_string.json')
expected_somatic_evidence_string = os.path.join(test_dir, 'resources/expected_somatic_evidence_string.json')


def get_test_clinvar_record():
    """The test file contains an extract of ClinVar XML for the record RCV000002127."""
    return [r for r in clinvar_xml_utils.ClinVarDataset(test_clinvar_record_file)][0]
