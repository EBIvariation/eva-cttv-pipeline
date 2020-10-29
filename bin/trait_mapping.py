import argparse
import sys

import eva_cttv_pipeline.trait_mapping.main as main


def launch():
    parser = ArgParser(sys.argv)

    main.main(parser.input_filepath, parser.output_mappings_filepath,
              parser.output_curation_filepath, parser.filters, parser.zooma_host,
              parser.oxo_target_list, parser.oxo_distance)


class ArgParser:

    def __init__(self, argv):
        description = """
                Script for running terms through Zooma, retrieving mapped uri, label from OLS,
                confidence of the mapping, and source of the mapping.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="input_filepath", required=True,
                            help="ClinVar XML dump file. One record per line.")
        parser.add_argument("-o", dest="output_mappings_filepath", required=True,
                            help="path to output file for mappings")
        parser.add_argument("-c", dest="output_curation_filepath", required=True,
                            help="path to output file for curation")
        parser.add_argument("-n", dest="ontologies", default="efo,ordo,hp,mondo",
                            help="ontologies to use in query")
        parser.add_argument("-r", dest="required", default="cttv,eva-clinvar,clinvar-xrefs,gwas",
                            help="data sources to use in query.")
        parser.add_argument("-p", dest="preferred", default="eva-clinvar,cttv,gwas,clinvar-xrefs",
                            help="preference for data sources, with preferred data source first.")
        parser.add_argument("-z", dest="zooma_host", default="https://www.ebi.ac.uk",
                            help="the host to use for querying zooma")
        parser.add_argument("-t", dest="oxo_target_list", default="Orphanet,efo,hp,mondo",
                            help="target ontologies to use with OxO")
        parser.add_argument("-d", dest="oxo_distance", default=3,
                            help="distance to use to query OxO.")

        args = parser.parse_args(args=argv[1:])

        self.input_filepath = args.input_filepath
        self.output_mappings_filepath = args.output_mappings_filepath
        self.output_curation_filepath = args.output_curation_filepath

        self.filters = {"ontologies": args.ontologies,
                        "required": args.required,
                        "preferred": args.preferred}

        self.zooma_host = args.zooma_host
        self.oxo_target_list = [target.strip() for target in args.oxo_target_list.split(",")]
        self.oxo_distance = args.oxo_distance


if __name__ == '__main__':
    launch()
