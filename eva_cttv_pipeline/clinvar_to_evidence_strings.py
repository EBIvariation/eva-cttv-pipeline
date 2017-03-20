import itertools
import json
import sys
from collections import defaultdict
from types import SimpleNamespace

import progressbar
import jsonschema

from eva_cttv_pipeline import cellbase_records
from eva_cttv_pipeline import efo_term
from eva_cttv_pipeline import config
from eva_cttv_pipeline import evidence_strings
from eva_cttv_pipeline import clinvar
from eva_cttv_pipeline import utilities
from eva_cttv_pipeline import consequence_type as CT
from eva_cttv_pipeline.trait import Trait


class Report:

    """
    Holds counters and other records of a pipeline run. Includes the list of evidence strings
    generated in the running of the pipeline.
    Includes method to write to output files, and __str__ shows the summary of the report.
    One instance of this class is instantiated in the running of the pipeline.
    """

    def __init__(self, unavailable_efo=None):
        if unavailable_efo is None:
            self.unavailable_efo = set()
        else:
            self.unavailable_efo = unavailable_efo

        self.unrecognised_clin_sigs = set()
        self.ensembl_gene_id_uris = set()
        self.traits = set()
        self.n_unrecognised_allele_origin = defaultdict(int)
        self.nsv_list = []
        self.unmapped_traits = defaultdict(int)
        self.evidence_string_list = []
        self.evidence_list = []  # To store Helen Parkinson records of the form
        self.counters = self.__get_counters()

    def __str__(self):

        report_strings = [
            str(self.counters["record_counter"]) + ' ClinVar records in total',
            str(len(self.evidence_string_list)) + ' evidence string jsons generated',
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
            str(self.counters["no_variant_to_ensg_mapping"]) +
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
            str(len(self.unavailable_efo)) +
            ' evidence strings were generated with traits without EFO correspondence',
            str(self.counters["n_valid_rs_and_nsv"]) +
            ' evidence strings were generated from ClinVar records with rs and nsv ids',
            str(self.counters["n_nsvs"]) + ' total nsvs found',
            str(self.counters["n_nsv_skipped_clin_sig"]) +
            ' ClinVar nsvs were skipped because of a different clinical significance',
            str(self.counters["n_nsv_skipped_wrong_ref_alt"]) +
            ' ClinVar nsvs were skipped because of same ref and alt'
        ])

        return '\n'.join(report_strings)

    def add_evidence_string(self, ev_string, clinvar_record, trait, ensembl_gene_id):
        try:
            ev_string.validate()
            self.evidence_string_list.append(ev_string)
        except jsonschema.exceptions.ValidationError as err:
            print('Error: evidence_string does not validate against schema.')
            # print('ClinVar accession: ' + record.clinvarRecord.accession)
            print(err)
            print(json.dumps(ev_string))
            print("clinvar record:\n%s" % clinvar_record)
            print("trait:\n%s" % trait)
            print("ensembl gene id: %s" % ensembl_gene_id)
            sys.exit(1)
        except efo_term.EFOTerm.IsObsoleteException as err:
            print('Error: obsolete EFO term.')
            print('Term: ' + ev_string.get_disease().efoid)
            print(err)
            print(json.dumps(ev_string))
            sys.exit(1)

    def write_output(self, dir_out):
        write_string_list_to_file(self.nsv_list, dir_out + '/' + config.NSV_LIST_FILE)

        # Contains traits without a mapping in Gary's xls
        with utilities.open_file(dir_out + '/' + config.UNMAPPED_TRAITS_FILE_NAME, 'wt') as fdw:
            fdw.write('Trait\tCount\n')
            for trait_list in self.unmapped_traits:
                fdw.write(str(trait_list) + '\t' +
                          str(self.unmapped_traits[trait_list]) + '\n')

        # Contains urls provided by Gary which are not yet included within EFO
        with utilities.open_file(dir_out + '/' + config.UNAVAILABLE_EFO_FILE_NAME, 'wt') as fdw:
            fdw.write('Trait\tCount\n')
            for clinvar_name in self.unavailable_efo:
                fdw.write(clinvar_name + "\n")

        with utilities.open_file(dir_out + '/' + config.EVIDENCE_STRINGS_FILE_NAME, 'wt') as fdw:
            for evidence_string in self.evidence_string_list:
                fdw.write(json.dumps(evidence_string) + '\n')

        with utilities.open_file(dir_out + '/' + config.EVIDENCE_RECORDS_FILE_NAME, 'wt') as fdw:
            for evidence_record in self.evidence_list:
                evidence_record_to_output = ['.' if ele is None else ele for ele in evidence_record]
                fdw.write('\t'.join(evidence_record_to_output) + '\n')

    @staticmethod
    def __get_counters():
        return {"n_processed_clinvar_records": 0,
                "n_pathogenic_no_rs": 0,
                "n_multiple_evidence_strings": 0,
                "n_multiple_allele_origin": 0,
                "n_germline_somatic": 0,
                "n_records_no_recognised_allele_origin": 0,
                "no_variant_to_ensg_mapping": 0,
                "n_more_than_one_efo_term": 0,
                "n_same_ref_alt": 0,
                "n_missed_strings_unmapped_traits": 0,
                "n_nsvs": 0,
                "n_valid_rs_and_nsv": 0,
                "n_nsv_skipped_clin_sig": 0,
                "n_nsv_skipped_wrong_ref_alt": 0,
                "record_counter": 0,
                "n_total_clinvar_records": 0}


def launch_pipeline(dir_out, allowed_clinical_significance, efo_mapping_file,
                    snp_2_gene_file, json_file):

    allowed_clinical_significance = allowed_clinical_significance.split(',') if \
        allowed_clinical_significance else get_default_allowed_clinical_significance()

    mappings = get_mappings(efo_mapping_file, snp_2_gene_file)

    report = clinvar_to_evidence_strings(allowed_clinical_significance, mappings, json_file)

    output(report, dir_out)


def output(report, dir_out):
    report.write_output(dir_out)
    print(report)


def file_len(fpath):
    i = -1
    with utilities.open_file(fpath, mode="rt") as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def clinvar_to_evidence_strings(allowed_clinical_significance, mappings, json_file):

    report = Report(mappings.unavailable_efo)

    cell_recs = cellbase_records.CellbaseRecords(json_file=json_file)

    num_cell_recs = file_len(json_file)

    bar = progressbar.ProgressBar(max_value=num_cell_recs)
    for cellbase_record in bar(cell_recs):
        n_ev_strings_per_record = 0
        clinvar_record = clinvar.ClinvarRecord(cellbase_record['clinvarSet'])

        for clinvar_record_measure in clinvar_record.measures:
            report.counters["record_counter"] += 1
            report.counters["n_nsvs"] += (clinvar_record_measure.nsv_id is not None)
            append_nsv(report.nsv_list, clinvar_record_measure)

            report.counters["n_multiple_allele_origin"] += (len(clinvar_record.allele_origins) > 1)

            traits = create_traits(clinvar_record.traits, mappings.trait_2_efo, report)

            converted_allele_origins = convert_allele_origins(clinvar_record.allele_origins)

            for consequence_type, trait, allele_origin in itertools.product(
                    get_consequence_types(clinvar_record_measure, mappings.consequence_type_dict),
                    traits,
                    converted_allele_origins):

                if skip_record(clinvar_record, clinvar_record_measure, consequence_type, allele_origin,
                               allowed_clinical_significance, report):
                    continue

                if allele_origin == 'germline':
                    evidence_string = evidence_strings.CTTVGeneticsEvidenceString(clinvar_record,
                                                                                  clinvar_record_measure,
                                                                                  report,
                                                                                  trait,
                                                                                  consequence_type)
                elif allele_origin == 'somatic':
                    evidence_string = evidence_strings.CTTVSomaticEvidenceString(clinvar_record,
                                                                                 clinvar_record_measure,
                                                                                 report,
                                                                                 trait,
                                                                                 consequence_type)
                report.add_evidence_string(evidence_string, clinvar_record, trait,
                                           consequence_type.ensembl_gene_id)
                report.evidence_list.append([clinvar_record.accession,
                                             clinvar_record_measure.rs_id,
                                             trait.clinvar_name,
                                             trait.ontology_id])
                report.counters["n_valid_rs_and_nsv"] += (clinvar_record_measure.nsv_id is not None)
                report.traits.add(trait.ontology_id)
                report.ensembl_gene_id_uris.add(
                    evidence_strings.get_ensembl_gene_id_uri(consequence_type.ensembl_gene_id))

                n_ev_strings_per_record += 1

            if n_ev_strings_per_record > 0:
                report.counters["n_processed_clinvar_records"] += 1
                if n_ev_strings_per_record > 1:
                    report.counters["n_multiple_evidence_strings"] += 1

    return report


def get_mappings(efo_mapping_file, snp_2_gene_file):
    mappings = SimpleNamespace()
    mappings.trait_2_efo, mappings.unavailable_efo = \
        load_efo_mapping(efo_mapping_file)

    mappings.consequence_type_dict = \
        CT.process_consequence_type_file(snp_2_gene_file)

    return mappings


def skip_record(clinvar_record, clinvar_record_measure, consequence_type, allele_origin,
                allowed_clinical_significance, report):
    if clinvar_record.clinical_significance not in allowed_clinical_significance:
        if clinvar_record_measure.nsv_id is not None:
            report.counters["n_nsv_skipped_clin_sig"] += 1
        return True

    if consequence_type is None:
        report.counters["no_variant_to_ensg_mapping"] += 1
        return True

    if allele_origin not in ('germline', 'somatic'):
        report.n_unrecognised_allele_origin[allele_origin] += 1
        return True

    return False


def get_consequence_types(clinvar_record_measure, consequence_type_dict):
    alt_str = clinvar_record_measure.alt if clinvar_record_measure.alt is not None else "-"
    coord_id = "{}:{}-{}:1/{}".format(clinvar_record_measure.chr, clinvar_record_measure.start,
                                      clinvar_record_measure.stop, alt_str)

    consequence_type_dict_id = None

    if clinvar_record_measure.rs_id is not None and \
                    clinvar_record_measure.rs_id in consequence_type_dict:
        consequence_type_dict_id = clinvar_record_measure.rs_id
    elif clinvar_record_measure.nsv_id is not None and \
                    clinvar_record_measure.nsv_id in consequence_type_dict:
        consequence_type_dict_id = clinvar_record_measure.nsv_id
    elif coord_id in consequence_type_dict:
        consequence_type_dict_id = coord_id
    elif clinvar_record_measure.clinvar_record.accession in consequence_type_dict:
        consequence_type_dict_id = clinvar_record_measure.clinvar_record.accession  # todo change this depending upon OT gene mapping file

    return consequence_type_dict[consequence_type_dict_id] \
        if consequence_type_dict_id is not None else [None]


def create_traits(clinvar_traits, trait_2_efo_dict, report):
    traits = []
    for trait_counter, name_list in enumerate(clinvar_traits):
        new_trait = create_trait(trait_counter, name_list, trait_2_efo_dict)
        if new_trait:
            traits.append(new_trait)
        else:
            report.counters["n_missed_strings_unmapped_traits"] += 1
            for name in name_list:
                report.unmapped_traits[name] += 1
    return traits


def create_trait(trait_counter, name_list, trait_2_efo_dict):
    trait = Trait(name_list, trait_counter, trait_2_efo_dict)
    # Only ClinVar records associated to a
    # trait with mapped EFO term will generate evidence_strings
    if trait.ontology_id is None:
        return None
    return trait


def write_string_list_to_file(string_list, filename):
    with utilities.open_file(filename, 'wt') as out_file:
        out_file.write('\n'.join(string_list))


def append_nsv(nsv_list, clinvar_record_measure):
    nsv = clinvar_record_measure.nsv_id
    if nsv is not None:
        nsv_list.append(nsv)
    return nsv_list


def load_efo_mapping(efo_mapping_file):
    trait_2_efo = {}
    unavailable_efo = set()
    n_efo_mappings = 0

    with utilities.open_file(efo_mapping_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line_list = line.rstrip().split("\t")
            clinvar_name = line_list[0].lower()
            if len(line_list) > 1:
                ontology_id = line_list[1]
                ontology_label = line_list[2] if len(line_list) > 2 else None
                trait_2_efo[clinvar_name] = (ontology_id, ontology_label)
                n_efo_mappings += 1
            else:
                unavailable_efo.add(clinvar_name)

    print(str(n_efo_mappings) + ' EFO mappings loaded')
    print(str(len(unavailable_efo)) + ' urls without an actual valid EFO mapping')

    return trait_2_efo, unavailable_efo


class UnhandledUrlTypeException(Exception):
    pass


def get_terms_from_file(terms_file_path):
    if terms_file_path is not None:
        print('Loading list of terms...')
        with utilities.open_file(terms_file_path, 'rt') as terms_file:
            terms_list = [line.rstrip() for line in terms_file]
        print(str(len(terms_file_path)) + ' terms found at ' + terms_file_path)
    else:
        terms_list = []

    return terms_list


def get_default_allowed_clinical_significance():
    return ['unknown', 'untested', 'non-pathogenic', 'probable-non-pathogenic',
            'probable-pathogenic', 'pathogenic', 'drug-response', 'drug response',
            'histocompatibility', 'other', 'benign', 'protective', 'not provided',
            'likely benign', 'confers sensitivity', 'uncertain significance',
            'likely pathogenic', 'conflicting data from submitters', 'risk factor',
            'association']


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



