import logging
import itertools
import json
import re
import sys
import os
from collections import defaultdict, Counter

import jsonschema

from eva_cttv_pipeline import clinvar_xml_utils
from eva_cttv_pipeline.evidence_string_generation import consequence_type as CT

logger = logging.getLogger(__package__)

# A regular expression to detect alleles with IUPAC ambiguity bases.
IUPAC_AMBIGUOUS_SEQUENCE = re.compile(r'[^ACGT]')

# Output file names.
EVIDENCE_STRINGS_FILE_NAME = 'evidence_strings.json'
EVIDENCE_RECORDS_FILE_NAME = 'evidence_records.tsv'
UNMAPPED_TRAITS_FILE_NAME = 'unmapped_traits.tsv'


class Report:
    """Holds counters and other records for a pipeline run."""

    def __init__(self, trait_mappings, consequence_mappings):
        self.report_strings = []

        # The main evidence string counter.
        self.evidence_string_count = 0

        # ClinVar record counters.
        self.clinvar_total = 0
        self.clinvar_fatal_no_allele_origin = 0
        self.clinvar_fatal_no_valid_traits = 0
        self.clinvar_skip_unsupported_variation = 0
        self.clinvar_skip_no_functional_consequences = 0
        self.clinvar_skip_missing_efo_mapping = 0
        self.clinvar_done_one_evidence_string = 0
        self.clinvar_done_multiple_evidence_strings = 0

        # Total number of trait-to-ontology mappings present in the database.
        self.total_trait_mappings = sum([len(mappings) for mappings in trait_mappings.values()])
        # All distinct (trait name, EFO ID) mappings used in the evidence strings.
        self.used_trait_mappings = set()
        # All unmapped trait names which prevented evidence string generation and their counts.
        self.unmapped_trait_names = Counter()

        # Variant-to-consequence mapping counts.
        self.total_consequence_mappings = sum([len(mappings) for mappings in consequence_mappings.values()])
        self.repeat_expansion_variants = 0

    def collate_report(self):
        # ClinVar tallies.
        clinvar_fatal = self.clinvar_fatal_no_allele_origin + self.clinvar_fatal_no_valid_traits
        clinvar_skipped = (self.clinvar_skip_unsupported_variation + self.clinvar_skip_no_functional_consequences +
                           self.clinvar_skip_missing_efo_mapping)
        clinvar_done = self.clinvar_done_one_evidence_string + self.clinvar_done_multiple_evidence_strings
        assert clinvar_fatal + clinvar_skipped + clinvar_done == self.clinvar_total, \
            'ClinVar evidence string tallies do not add up to the total amount.'

        return f'''Total number of evidence strings generated\t{self.evidence_string_count}

            Total number of ClinVar records\t{self.clinvar_total}
                Fatal: Cannot be processed ever\t{clinvar_fatal}
                    No allele origin\t{self.clinvar_fatal_no_allele_origin}
                    No traits with valid names\t{self.clinvar_fatal_no_valid_traits}
                Skipped: Can be rescued by future improvements\t{clinvar_skipped}
                    Unsupported variation type\t{self.clinvar_skip_unsupported_variation}
                    No functional consequences\t{self.clinvar_skip_no_functional_consequences}
                    Missing EFO mapping\t{self.clinvar_skip_missing_efo_mapping}
                Done: Generated at least one evidence string\t{clinvar_done}
                    One evidence string\t{self.clinvar_done_one_evidence_string}
                    Multiple evidence strings\t{self.clinvar_done_multiple_evidence_strings}
            Percentage of all potentially supportable ClinVar records which generated at least one evidence string\t{
                clinvar_done / (clinvar_skipped + clinvar_done):.1%}

            Total number of trait-to-ontology mappings in the database\t{self.total_trait_mappings}
                The number of distinct trait-to-ontology mappings used in the evidence strings\t{
                    len(self.used_trait_mappings)}
            The number of distinct unmapped trait names which prevented evidence string generation\t{
                len(self.unmapped_trait_names)}

            Total number of variant to consequence mappings\t{self.total_consequence_mappings}
                Number of repeat expansion variants\t{self.repeat_expansion_variants}'''.replace('\n' + ' ' * 12, '\n')

    def write_unmapped_terms(self, dir_out):
        with open(os.path.join(dir_out, UNMAPPED_TRAITS_FILE_NAME), 'w') as unmapped_traits_file:
            for trait, number_of_occurrences in sorted(self.unmapped_trait_names.items(), key=lambda x: -x[1]):
                unmapped_traits_file.write(f'{trait}\t{number_of_occurrences}\n')


def validate_evidence_string(ev_string, ot_schema_contents):
    try:
        jsonschema.validate(ev_string, ot_schema_contents, format_checker=jsonschema.FormatChecker())
    except jsonschema.exceptions.ValidationError as err:
        logger.error('Error: evidence string does not validate against schema.')
        logger.error(f'Error message: {err}')
        logger.error(f'Complete evidence string: {json.dumps(ev_string)}')
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
        output_evidence_strings=os.path.join(dir_out, EVIDENCE_STRINGS_FILE_NAME))
    print(report.collate_report())
    report.write_unmapped_terms(dir_out)


def clinvar_to_evidence_strings(string_to_efo_mappings, variant_to_gene_mappings, clinvar_xml, ot_schema,
                                output_evidence_strings):
    report = Report(trait_mappings=string_to_efo_mappings, consequence_mappings=variant_to_gene_mappings)
    ot_schema_contents = json.loads(open(ot_schema).read())
    output_evidence_strings_file = open(output_evidence_strings, 'wt')

    logger.info('Processing ClinVar records')
    for clinvar_record in clinvar_xml_utils.ClinVarDataset(clinvar_xml):
        report.clinvar_total += 1
        if report.clinvar_total % 1000 == 0:
            logger.info(f'{report.clinvar_total} records processed')
        if clinvar_record.measure and clinvar_record.measure.is_repeat_expansion_variant:
            report.repeat_expansion_variants += len(get_consequence_types(clinvar_record.measure,
                                                                          variant_to_gene_mappings))

        # Failure mode 1 (fatal). A ClinVar record contains no allele origin.
        if not clinvar_record.allele_origins:
            report.clinvar_fatal_no_allele_origin += 1
            continue

        # Failure mode 2 (fatal). A ClinVar record contains no valid traits (traits which have at least one valid,
        # potentially mappable name).
        if not clinvar_record.traits_with_valid_names:
            report.clinvar_fatal_no_valid_traits += 1
            continue

        # Failure mode 3 (skip). A ClinVar record contains an unsupported variation type.
        if clinvar_record.measure is None:
            report.clinvar_skip_unsupported_variation += 1
            continue

        # Within each ClinVar record, an evidence string is generated for all possible permutations of (1) allele
        # origins, (2) EFO mappings, and (3) genes where the variant has effect.
        grouped_allele_origins = convert_allele_origins(clinvar_record.allele_origins)
        consequence_types = get_consequence_types(clinvar_record.measure, variant_to_gene_mappings)
        grouped_diseases = group_diseases_by_efo_mapping(clinvar_record.traits_with_valid_names,
                                                         string_to_efo_mappings)

        # Failure mode 4 (skip). No functional consequences are available.
        if not consequence_types:
            report.clinvar_skip_no_functional_consequences += 1
            continue

        # Failure mode 5 (skip). A ClinVar records has at least one trait with at least one valid name, but no suitable
        # EFO mappings were found in the database.
        if not grouped_diseases:
            report.clinvar_skip_missing_efo_mapping += 1
            unmapped_trait_name = clinvar_record.traits_with_valid_names[0].preferred_or_other_valid_name
            report.unmapped_trait_names[unmapped_trait_name] += 1
            continue

        assert grouped_allele_origins and grouped_diseases and consequence_types, \
            'Some of the attribute lists are still empty even after passing all checks.'

        evidence_strings_generated = 0
        for allele_origins, disease_attributes, consequence_attributes in itertools.product(
                grouped_allele_origins, grouped_diseases, consequence_types):
            disease_name, disease_source_id, disease_mapped_efo_id = disease_attributes
            evidence_string = generate_evidence_string(clinvar_record, allele_origins, disease_name, disease_source_id,
                                                       disease_mapped_efo_id, consequence_attributes)

            # Validate and immediately output the evidence string (not keeping everything in memory).
            validate_evidence_string(evidence_string, ot_schema_contents)
            output_evidence_strings_file.write(json.dumps(evidence_string) + '\n')

            # Record some evidence string and trait metrics.
            evidence_strings_generated += 1
            report.used_trait_mappings.add((disease_name, disease_mapped_efo_id))

        assert evidence_strings_generated != 0, 'No evidence strings generated despite all attributes passing checks.'
        if evidence_strings_generated == 1:
            report.clinvar_done_one_evidence_string += 1
        else:
            report.clinvar_done_multiple_evidence_strings += 1
        report.evidence_string_count += evidence_strings_generated

    output_evidence_strings_file.close()
    return report


def generate_evidence_string(clinvar_record, allele_origins, disease_name, disease_source_id, disease_mapped_efo_id,
                             consequence_attributes):
    """Generates an evidence string based on ClinVar record and some additional attributes."""
    is_somatic = allele_origins == ['somatic']
    evidence_string = {
        # ALLELE ORIGIN ATTRIBUTES. There are three attributes which are currently completely redundant (their
        # values are completely correlated between each other). This is intended, see answer to question 1 here:
        # https://github.com/EBIvariation/eva-opentargets/issues/189#issuecomment-782136128.
        'alleleOrigins': allele_origins,
        'datasourceId': 'eva_somatic' if is_somatic else 'eva',
        'datatypeId': 'somatic_mutation' if is_somatic else 'genetic_association',

        # ASSOCIATION ATTRIBUTES.
        # List of patterns of inheritance reported for the variant.
        'allelicRequirements': clinvar_record.mode_of_inheritance,

        # Levels of clinical significance reported for the variant.
        'clinicalSignificances': clinvar_record.clinical_significance_list,

        # Confidence (review status).
        'confidence': clinvar_record.review_status,

        # Literature. ClinVar records provide three types of references: trait-specific; variant-specific; and
        # "observed in" references. Open Targets are interested only in that last category.
        'literature': sorted(set([str(r) for r in clinvar_record.evidence_support_pubmed_refs])),

        # RCV identifier.
        'studyId': clinvar_record.accession,

        # VARIANT ATTRIBUTES.
        'targetFromSourceId': consequence_attributes.ensembl_gene_id,
        'variantFunctionalConsequenceId': consequence_attributes.so_term.accession,
        'variantId': clinvar_record.measure.vcf_full_coords,  # CHROM_POS_REF_ALT notation.
        'variantRsId': clinvar_record.measure.rs_id,

        # PHENOTYPE ATTRIBUTES.
        # The alphabetical list of *all* disease names from that ClinVar record.
        'cohortPhenotypes': sorted([trait.preferred_or_other_valid_name for trait in clinvar_record.traits
                                    if trait.preferred_or_other_valid_name is not None]),

        # One disease name for this evidence string (see group_diseases_by_efo_mapping).
        'diseaseFromSource': disease_name,

        # The internal identifier of that disease.
        'diseaseFromSourceId': disease_source_id,

        # The EFO identifier to which we mapped that first disease. Converting the URI to a compact representation as
        # required by the Open Targets JSON schema.
        'diseaseFromSourceMappedId': disease_mapped_efo_id.split('/')[-1],
    }
    # Remove the attributes with empty values (either None or empty lists).
    evidence_string = {key: value for key, value in evidence_string.items() if value}
    return evidence_string


def get_consequence_types(clinvar_record_measure, consequence_type_dict):
    """Returns the list of functional consequences for a given ClinVar record measure.

    This is the place where ClinVar records are paired with the information about gene and functional consequences.
    This information is produced by two pipelines in the `consequence_prediction` subdirectory:
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
        # Example of such an identifier: 14:23423715:G:A.
        coord_id = ':'.join([clinvar_record_measure.chr, str(clinvar_record_measure.vcf_pos),
                             clinvar_record_measure.vcf_ref, clinvar_record_measure.vcf_alt])
        # Log unusual variants with IUPAC ambiguity symbols. Example: 12_32625716_G_H (from RCV000032000).
        if IUPAC_AMBIGUOUS_SEQUENCE.search(clinvar_record_measure.vcf_ref + clinvar_record_measure.vcf_alt):
            logger.warning(f'Observed variant with non-ACGT allele sequences: {coord_id}')
        if coord_id in consequence_type_dict:
            return consequence_type_dict[coord_id]

    # Previously, the pairing was also attempted based on rsID and nsvID. This is not reliable because of lack of allele
    # specificity, and has been removed.
    return []


def write_string_list_to_file(string_list, filename):
    with open(filename, 'wt') as out_file:
        out_file.write('\n'.join(string_list))


def load_efo_mapping(efo_mapping_file):
    trait_2_efo = defaultdict(list)
    n_efo_mappings = 0

    with open(efo_mapping_file, 'rt') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#') or not line:
                continue
            line_list = line.split('\t')
            assert len(line_list) == 3, f'Incorrect string to EFO mapping format for line {line}'
            clinvar_name, ontology_id, ontology_label = line_list
            trait_2_efo[clinvar_name.lower()].append((ontology_id, ontology_label))
            n_efo_mappings += 1
    logger.info('{} EFO mappings loaded'.format(n_efo_mappings))
    return trait_2_efo


def get_terms_from_file(terms_file_path):
    if terms_file_path is not None:
        print('Loading list of terms...')
        with open(terms_file_path, 'rt') as terms_file:
            terms_list = [line.rstrip() for line in terms_file]
        print(str(len(terms_file_path)) + ' terms found at ' + terms_file_path)
    else:
        terms_list = []

    return terms_list


def convert_allele_origins(orig_allele_origins):
    """Splits the original list of allele origins from ClinVar into up to two groups: one for 'somatic' (if present),
    and another one for all other types, which are considered 'germline' (if present)."""
    orig_allele_origins = {item.lower() for item in orig_allele_origins}
    converted_allele_origins = []
    if 'somatic' in orig_allele_origins:
        converted_allele_origins.append(['somatic'])
        orig_allele_origins.remove('somatic')
    if orig_allele_origins:  # Is there something remaining for the second (non-somatic) group?
        converted_allele_origins.append(sorted(orig_allele_origins))
    return converted_allele_origins


def group_diseases_by_efo_mapping(clinvar_record_traits, string_to_efo_mappings):
    """Processes and groups diseases from a ClinVar record in two ways. First, diseases mapping to the same EFO term are
    grouped together. Second, if a group of diseases has multiple EFO mappings, they are split into separate groups.
    For each group, the tuple of the following values is returned: (1) the original name of the disease which comes
    first lexicographically; (2) the original MedGen ID of that disease; (3) the mapped EFO term of that disease.
    Diseases without mappings are skipped.

    For example, for a ClinVar record with the following mappings:
        * Diseases A, B, C -> EFO_1
        * Disease D -> EFO_2 & EFO_3
        * Diseases E, F -> EFO_4 & EFO_5
        * Disease G -> No mapping

    The following tuples will be generated:
        * (A, MedGen_A, EFO_1)
        * (D, MedGen_D, EFO_2)
        * (D, MedGen_D, EFO_3)
        * (E, MedGen_E, EFO_4)
        * (E, MedGen_E, EFO_5)"""

    # Group traits by their EFO mappings and explode multiple mappings.
    efo_to_traits = defaultdict(list)  # Key: EFO ID, value: list of traits mapped to that ID.
    for trait in clinvar_record_traits:
        # Try to match using all trait names.
        for trait_name in trait.all_names:
            for efo_id, efo_label in string_to_efo_mappings.get(trait_name.lower(), []):
                efo_to_traits[efo_id].append(trait)

    # Generate tuples by keeping only one disease from each group.
    grouped_tuples = []
    for efo_id, traits in efo_to_traits.items():
        traits = sorted(traits, key=lambda t: t.preferred_or_other_valid_name)
        selected_trait = traits[0]
        grouped_tuples.append((selected_trait.preferred_or_other_valid_name, selected_trait.medgen_id, efo_id))
    return grouped_tuples
