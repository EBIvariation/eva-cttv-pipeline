import os

from eva_cttv_pipeline import clinvar_xml_utils

test_dir = os.path.dirname(__file__)
efo_mapping_file = os.path.join(test_dir, 'resources', 'feb16_jul16_combined_trait_to_url.tsv')
snp_2_gene_file = os.path.join(test_dir, 'resources/snp2gene_extract.tsv')
test_clinvar_record_file = os.path.join(test_dir, 'resources/test_clinvar_record.xml.gz')
open_targets_schema_gz = os.path.join(test_dir, 'resources/opentargets.1.7.3.json.gz')
expected_genetics_evidence_string = os.path.join(test_dir, 'resources/expected_genetics_evidence_string.json')
expected_somatic_evidence_string = os.path.join(test_dir, 'resources/expected_somatic_evidence_string.json')


def get_test_clinvar_record():
    """The test file contains an extract of ClinVar XML for the record RCV000002127."""
    return [r for r in clinvar_xml_utils.ClinVarDataset(test_clinvar_record_file)][0]
