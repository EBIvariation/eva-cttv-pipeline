#!/usr/bin/python

import sys

from eva_cttv_pipeline import utilities, clinvar_to_evidence_strings


def main():
    parser = utilities.ArgParser(sys.argv)

    utilities.check_for_local_schema()

    utilities.check_dir_exists_create(parser.out)

    clinvar_to_evidence_strings.launch_pipeline(parser.out,
                                                allowed_clinical_significance=parser.clinical_significance,
                                                ignore_terms_file=parser.ignore_terms_file,
                                                adapt_terms_file=parser.adapt_terms_file,
                                                efo_mapping_file=parser.efo_mapping_file,
                                                snp_2_gene_file=parser.snp_2_gene_file,
                                                variant_summary_file=parser.variant_summary_file,
                                                json_file=parser.json_file)

    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Finished <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')


if __name__ == '__main__':
    main()
