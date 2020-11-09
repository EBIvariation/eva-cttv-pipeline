import logging
import itertools
import copy
import json
import sys
import os
from collections import defaultdict

import jsonschema

from eva_cttv_pipeline import clinvar_xml_utils, file_utils
from eva_cttv_pipeline.evidence_string_generation import config
from eva_cttv_pipeline.evidence_string_generation import evidence_strings
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT

logger = logging.getLogger(__package__)


class Report:
    """Holds counters and other records of a pipeline run. Includes method to write to output files, and __str__ shows
    the summary of the report. One instance of this class is instantiated in the running of the pipeline."""

    def __init__(self, trait_mappings=None):
        if trait_mappings is None:
            self.trait_mappings = {}
        else:
            self.trait_mappings = copy.deepcopy(trait_mappings)

        self.unrecognised_clin_sigs = set()
        self.ensembl_gene_id_uris = set()
        self.traits = set()
        self.n_unrecognised_allele_origin = defaultdict(int)
        self.nsv_list = []
        self.unmapped_traits = defaultdict(int)
        self.evidence_string_count = 0
        self.counters = self.__get_counters()

    def __str__(self):

        report_strings = [
            str(self.counters["record_counter"]) + ' ClinVar records in total',
            str(self.evidence_string_count) + ' evidence string jsons generated',
            str(self.counters["n_processed_clinvar_records"]) +
            ' ClinVar records generated at least one evidence string',
            str(len(self.unrecognised_clin_sigs)) +
            " Clinical significance string(s) not found " +
            "among those described in ClinVar documentation:",
            str(self.unrecognised_clin_sigs),
            str(self.counters["n_same_ref_alt"]) +
            ' ClinVar records with allowed clinical significance ' +
            'did present the same reference and alternate and were skipped',
            'Activities of those ClinVar records with ' +
            'unrecognized clinical significances were set to "unknown".',
            str(len(self.ensembl_gene_id_uris)) +
            ' distinct ensembl gene ids appear in generated evidence string json objects',
            str(len(self.traits)) +
            ' distinct trait names found to include in generated evidence string json objects',
            str(self.counters["n_pathogenic_no_rs"]) +
            ' ClinVar records with allowed clinical significance DO NOT have an rs id',
            str(self.counters["n_multiple_evidence_strings"]) +
            ' ClinVar records generated more than one evidence_string',
            str(self.counters["n_germline_somatic"]) +
            ' ClinVar records with germline and somatic origins',
            str(self.counters["n_multiple_allele_origin"]) +
            ' ClinVar records with more than one allele origin',
            'Number valid ClinVar records with unprocessed allele origins:'
        ]

        report_strings.extend(
            [' ' + alleleOrigin + ': ' +
             str(self.n_unrecognised_allele_origin[alleleOrigin])
             for alleleOrigin in self.n_unrecognised_allele_origin])

        report_strings.extend([
            str(self.counters["n_no_variant_to_ensg_mapping"]) +
            ' ClinVar records with allowed clinical significance and valid rs id ' +
            'were skipped due to a lack of Variant->ENSG mapping.',
            str(self.counters["n_missed_strings_unmapped_traits"]) +
            ' ClinVar records with allowed clinical significance, valid rs id and ' +
            'Variant->ENSG mapping were skipped due to a lack of EFO mapping (see ' +
            config.UNMAPPED_TRAITS_FILE_NAME + ').',
            str(self.counters["n_records_no_recognised_allele_origin"]) +
            ' ClinVar records with allowed clinical significance, ' +
            'valid rs id, valid Variant->ENSG' +
            ' mapping and valid EFO mapping were skipped due to a lack of a valid alleleOrigin.',
            str(self.counters["n_more_than_one_efo_term"]) +
            ' evidence strings with more than one trait mapped to EFO terms',
            str(self.counters["n_valid_rs_and_nsv"]) +
            ' evidence strings were generated from ClinVar records with rs and nsv ids',
            str(self.counters["n_nsvs"]) + ' total nsvs found',
            str(self.counters["n_nsv_skipped_clin_sig"]) +
            ' ClinVar nsvs were skipped because of a different clinical significance',
            str(self.counters["n_nsv_skipped_wrong_ref_alt"]) +
            ' ClinVar nsvs were skipped because of same ref and alt'
        ])

        return '\n'.join(report_strings)

    def write_output(self, dir_out):
        write_string_list_to_file(self.nsv_list, dir_out + '/' + config.NSV_LIST_FILE)

        # Contains traits without a mapping in Gary's xls
        with file_utils.open_file(dir_out + '/' + config.UNMAPPED_TRAITS_FILE_NAME, 'wt') as fdw:
            fdw.write('Trait\tCount\n')
            for trait_list in self.unmapped_traits:
                fdw.write(str(trait_list) + '\t' +
                          str(self.unmapped_traits[trait_list]) + '\n')

    def remove_trait_mapping(self, trait_name):
        if trait_name in self.trait_mappings:
            del self.trait_mappings[trait_name]

    @staticmethod
    def __get_counters():
        return {"n_processed_clinvar_records": 0,
                "n_pathogenic_no_rs": 0,
                "n_multiple_evidence_strings": 0,
                "n_multiple_allele_origin": 0,
                "n_germline_somatic": 0,
                "n_records_no_recognised_allele_origin": 0,
                "n_no_variant_to_ensg_mapping": 0,
                "n_more_than_one_efo_term": 0,
                "n_same_ref_alt": 0,
                "n_missed_strings_unmapped_traits": 0,
                "n_nsvs": 0,
                "n_valid_rs_and_nsv": 0,
                "n_nsv_skipped_clin_sig": 0,
                "n_nsv_skipped_wrong_ref_alt": 0,
                "record_counter": 0,
                "n_total_clinvar_records": 0}


def validate_evidence_string(ev_string, clinvar_record, trait, ensembl_gene_id, ot_schema_contents):
    try:
        ev_string.validate(ot_schema_contents)
    except jsonschema.exceptions.ValidationError as err:
        logger.error('Error: evidence string does not validate against schema.')
        logger.error(f'Error message: {err}')
        logger.error(f'Complete evidence string: {json.dumps(ev_string)}')
        logger.error(f'ClinVar record: {clinvar_record}')
        logger.error(f'ClinVar trait: {trait}')
        logger.error(f'Ensembl gene ID: {ensembl_gene_id}')
        sys.exit(1)
    except jsonschema.exceptions.SchemaError:
        logger.error('Error: OpenTargets schema file is invalid')
        sys.exit(1)


def launch_pipeline(clinvar_xml_file, efo_mapping_file, gene_mapping_file, ot_schema_file, dir_out):
    os.makedirs(dir_out, exist_ok=True)
    string_to_efo_mappings = load_efo_mapping(efo_mapping_file)
    variant_to_gene_mappings = CT.process_consequence_type_file(gene_mapping_file)
    report = clinvar_to_evidence_strings(
        string_to_efo_mappings, variant_to_gene_mappings, clinvar_xml_file, ot_schema_file,
        output_evidence_strings=os.path.join(dir_out, config.EVIDENCE_STRINGS_FILE_NAME))
    report.write_output(dir_out)
    print(report)


def clinvar_to_evidence_strings(string_to_efo_mappings, variant_to_gene_mappings, clinvar_xml, ot_schema,
                                output_evidence_strings):
    report = Report(trait_mappings=string_to_efo_mappings)
    ot_schema_contents = json.loads(open(ot_schema).read())
    output_evidence_strings_file = file_utils.open_file(output_evidence_strings, 'wt')

    logger.info('Processing ClinVar records')
    for clinvar_record in clinvar_xml_utils.ClinVarDataset(clinvar_xml):
        report.counters['record_counter'] += 1
        if report.counters['record_counter'] % 1000 == 0:
            logger.info('{} records processed'.format(report.counters['record_counter']))
        n_ev_strings_per_record = 0

        if clinvar_record.measure is None:  # No valid variants to process in this record
            continue
        report.counters['n_nsvs'] += (clinvar_record.measure.nsv_id is not None)
        append_nsv(report.nsv_list, clinvar_record.measure)
        report.counters['n_multiple_allele_origin'] += (len(clinvar_record.allele_origins) > 1)
        converted_allele_origins = convert_allele_origins(clinvar_record.allele_origins)
        for consequence_type, clinvar_trait, allele_origin in itertools.product(
                get_consequence_types(clinvar_record.measure, variant_to_gene_mappings),
                clinvar_record.traits,
                converted_allele_origins):

            # Skip records for which crucial information (consequences or EFO mappings) is not available
            if consequence_type is None:
                report.counters['n_no_variant_to_ensg_mapping'] += 1
                continue
            if clinvar_trait.name.lower() not in string_to_efo_mappings:
                report.counters['n_missed_strings_unmapped_traits'] += 1
                continue

            # Iterate over all EFO mappings (there may be multiple per trait)
            for ontology_id, ontology_label in string_to_efo_mappings[clinvar_trait.name.lower()]:
                if allele_origin == 'germline':
                    evidence_string = evidence_strings.CTTVGeneticsEvidenceString(
                        clinvar_record, clinvar_trait, ontology_id, ontology_label, consequence_type)
                elif allele_origin == 'somatic':
                    evidence_string = evidence_strings.CTTVSomaticEvidenceString(
                        clinvar_record, clinvar_trait, ontology_id, ontology_label, consequence_type)
                else:
                    raise AssertionError('Unknown allele_origin present in the data: {}'.format(allele_origin))

                # Validate and immediately output the evidence string (not keeping everything in memory)
                validate_evidence_string(evidence_string, clinvar_record, clinvar_trait,
                                         consequence_type.ensembl_gene_id, ot_schema_contents)
                output_evidence_strings_file.write(json.dumps(evidence_string) + '\n')
                report.evidence_string_count += 1
                report.counters['n_valid_rs_and_nsv'] += (clinvar_record.measure.nsv_id is not None)
                report.traits.add(ontology_id)
                report.remove_trait_mapping(clinvar_trait.name)
                report.ensembl_gene_id_uris.add(
                    evidence_strings.get_ensembl_gene_id_uri(consequence_type.ensembl_gene_id))
                n_ev_strings_per_record += 1

        if n_ev_strings_per_record > 0:
            report.counters['n_processed_clinvar_records'] += 1
            if n_ev_strings_per_record > 1:
                report.counters['n_multiple_evidence_strings'] += 1

    output_evidence_strings_file.close()
    return report


def get_consequence_types(clinvar_record_measure, consequence_type_dict):
    """Returns the list of functional consequences for a given ClinVar record measure.

    This is the place where ClinVar records are paired with the information about gene and functional consequences.
    This information is produced by two pipelines in the `vep-mapping-pipeline` subdirectory:
    1. The main one, `vep_mapping_pipeline`, runs all records in the ClinVar VCF dump through Variant Effect Predictor
       and outputs the results using a VCF-compatible "CHROM:POS:REF:ALT" identifier.
    2. The auxiliary one, `repeat_expansion_variants`, uses a different approach to extract information about repeat
       expansion variants from the ClinVar TSV dump. The repeat expansion variants have several important features:
       a. Most (but not all) of them are not present in ClinVar VCF dump, only in the TSV dump.
       b. Even the variants which have the coordinates cannot be adequately processed by VEP: it will output the wrong
          functional consequence type.
       c. The "CHROM:POS:REF:ALT" notation is not useful for these variants, because some of them have an indeterminate
          number of repeats, hence there is no single ALT allele.
       This second pipeline outputs the results using the RCV identifier from ClinVar."""

    # As the first option, try to link a variant and its functional consequence by RCV accession. This will only happen
    # for repeat expansion variants, since the main pipeline uses CHROM:POS:REF:ALT identifiers.
    #
    # The reason we don't use RCV accessions more often is because they are, in general, not variant-specific: the same
    # RCV may link to multiple variants with different functional consequences. However, for repeat expansion variants
    # RCV accessions *are* specific and contain only one variant.
    #
    # The reason we pair first by RCV accession and *then* by CHROM:POS:REF:ALT identifiers is that some repeat
    # expansion variants actually *do* appear in the ClinVar VCF dump, will be fed to VEP, and will receive incorrect
    # consequence annotations. By using RCV pairing first, we prioritise results of the variant expansion pipeline over
    # the general VEP pipeline.
    if clinvar_record_measure.clinvar_record.accession in consequence_type_dict:
        return consequence_type_dict[clinvar_record_measure.clinvar_record.accession]

    # If RCV is not present in the consequences file, pair using full variant description (CHROM:POS:REF:ALT)
    if clinvar_record_measure.has_complete_coordinates:
        # This VCF-flavoured identifier is used to pair ClinVar records with functional consequence predictions.
        # Example of such an identifier: 14:23423715:G:A
        coord_id = ':'.join([clinvar_record_measure.chr, str(clinvar_record_measure.vcf_pos),
                             clinvar_record_measure.vcf_ref, clinvar_record_measure.vcf_alt])
        if coord_id in consequence_type_dict:
            return consequence_type_dict[coord_id]

    # Previously, the pairing was also attempted based on rsID and nsvID. This is not reliable because of lack of allele
    # specificity, and has been removed.
    return [None]


def write_string_list_to_file(string_list, filename):
    with file_utils.open_file(filename, 'wt') as out_file:
        out_file.write('\n'.join(string_list))


def append_nsv(nsv_list, clinvar_record_measure):
    nsv = clinvar_record_measure.nsv_id
    if nsv is not None:
        nsv_list.append(nsv)
    return nsv_list


def load_efo_mapping(efo_mapping_file):
    trait_2_efo = defaultdict(list)
    n_efo_mappings = 0

    with file_utils.open_file(efo_mapping_file, "rt") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#") or not line:
                continue
            line_list = line.split("\t")
            clinvar_name = line_list[0].lower()
            if len(line_list) > 1:
                ontology_id_list = line_list[1].split("|")
                ontology_label_list = line_list[2].split("|") if len(line_list) > 2 else [None] * len(ontology_id_list)
                for ontology_id, ontology_label in zip(ontology_id_list, ontology_label_list):
                    trait_2_efo[clinvar_name].append((ontology_id, ontology_label))
                n_efo_mappings += 1
            else:
                raise ValueError('No mapping provided for trait: {}'.format(clinvar_name))
    logger.info('{} EFO mappings loaded'.format(n_efo_mappings))
    return trait_2_efo


def get_terms_from_file(terms_file_path):
    if terms_file_path is not None:
        print('Loading list of terms...')
        with file_utils.open_file(terms_file_path, 'rt') as terms_file:
            terms_list = [line.rstrip() for line in terms_file]
        print(str(len(terms_file_path)) + ' terms found at ' + terms_file_path)
    else:
        terms_list = []

    return terms_list


def convert_allele_origins(orig_allele_origins):
    orig_allele_origins = [item.lower() for item in orig_allele_origins]
    converted_allele_origins = []
    if "somatic" in orig_allele_origins:
        converted_allele_origins.append("somatic")
    if set(orig_allele_origins).intersection({"biparental", "de novo", "germline", "inherited",
                                              "maternal", "not applicable", "not provided",
                                              "paternal", "uniparental", "unknown"}):
        converted_allele_origins.append("germline")

    return converted_allele_origins
