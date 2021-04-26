#!/usr/bin/env python3
"""Pipeline for mapping variants to the genes they affect and their functional consequences, using Ensembl VEP API. For
documentation, refer to /README.md"""

import argparse
import itertools
import json
import logging
import os
import requests
import sys

from retry import retry

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    '--enable-distant-querying', action='store_true',
    help='Enables a second iteration of querying VEP for distant gene variants, which is disabled by default'
)
parser.add_argument(
    '--report-distance', action='store_true',
    help='Report a distance to the gene for upstream/downstream gene variants. Disabled by default'
)

logging.basicConfig()
logger = logging.getLogger('consequence_mapping')
logger.setLevel(logging.INFO)

# The "distance to the nearest gene" parameters, used to query VEP in first and second iterations, respectively.
VEP_SHORT_QUERY_DISTANCE = 5000
VEP_LONG_QUERY_DISTANCE = 500000


def deduplicate_list(lst):
    """Removes duplicates from a list containing arbitrary (possibly unhashable) values."""
    return [element for element, _ in itertools.groupby(sorted(lst))]


def colon_based_id_to_vep_id(colon_id):
    """Converts a colon-based identifier to VEP compatible one. Example: '15:7237571:C:T' → '15 7237571 . C T'"""
    id_fields = colon_id.split(':')
    assert len(id_fields) == 4, 'Invalid colon-based identifier supplied (should contain exactly 4 fields)'
    return '{} {} . {} {}'.format(*id_fields)


def vep_id_to_colon_id(vep_id):
    """Converts a specific type of VEP compatible identifier to colon-based one. VEP supports several types of variant
    identifiers. This function only converts a single type, of the form 'CHROM POS . REF ALT', delimited by spaces and
    compatible with the first five columns of VCF."""
    vep_id_fields = vep_id.split(' ')
    return ':'.join([vep_id_fields[0], vep_id_fields[1], vep_id_fields[3], vep_id_fields[4]])


@retry(tries=10, delay=5, backoff=1.2, jitter=(1, 3), logger=logger)
def query_vep(variants, search_distance):
    """Query VEP and return results in JSON format. Upstream/downstream genes are searched up to a given distance."""
    ensembl_request_url = 'https://rest.ensembl.org/vep/human/region'
    headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
    result = requests.post(ensembl_request_url, headers=headers, data=json.dumps({
        'variants': variants, 'distance': search_distance, 'shift_3prime': 0,
    }))
    # If there was an HTTP error, raise an exception. This will be caught by @retry
    result.raise_for_status()
    return result.json()


def load_consequence_severity_rank():
    """Loads severity rankings for consequence terms."""
    severity_ranking_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'severity_ranking.txt')
    severity_ranking = open(severity_ranking_path).read().splitlines()
    return {term: index for index, term in enumerate(severity_ranking)}


def extract_consequences(vep_results, acceptable_biotypes, only_closest, results_by_variant, report_distance=False):
    """Given VEP results, return a list of consequences matching certain criteria.

    Args:
        vep_results: results obtained from VEP for a list of variants, in JSON format.
        acceptable_biotypes: a list of transcript biotypes to consider (as defined in Ensembl documentation, see
            https://www.ensembl.org/info/genome/genebuild/biotypes.html). Consequences for other transcript biotypes
            will be ignored.
        only_closest: if this flag is specified, then at most one consequence per variant will be returned. The
            consequences are sorted by distance from the gene and the closest one is chosen. In case of a tie,
            consequence is selected at random. If this flag is not specified, all consequences for each variant will
            be returned.
        results_by_variant: a dict to write the results into.
        report_distance: if set to True, a distance from the variant to the gene it affects will be reported
            (applicable to upstream and downstream gene variants). Otherwise, the distance parameter will always be 0.
    """
    consequence_term_severity_rank = load_consequence_severity_rank()
    for result in vep_results:
        variant_identifier = result['input']
        results_by_variant.setdefault(variant_identifier, [])
        consequences = result.get('transcript_consequences', [])

        # Keep only consequences with the allowed biotypes; if there are none, skip this variant
        consequences = [c for c in consequences if c['biotype'] in acceptable_biotypes]
        if not consequences:
            continue

        # Flatten the list of consequence terms and find the most severe one
        all_consequence_terms = [term for c in consequences for term in c['consequence_terms']]
        all_consequence_terms.sort(key=lambda term: consequence_term_severity_rank[term])
        most_severe_consequence_term = all_consequence_terms[0]

        # Keep only consequences which include the most severe consequence term; sort by increasing order of distance.
        # If there is no 'distance' attribute in VEP results, it means that it is not applicable as the variant resides
        # *inside* the gene; hence, in this case the distance is set to 0.
        consequences = [c for c in consequences if most_severe_consequence_term in c['consequence_terms']]
        consequences.sort(key=lambda consequence: abs(consequence.get('distance', 0)))

        # If mandated by a flag, keep only one (least distant) consequence
        if only_closest:
            consequences = [consequences[0]]

        # Return a subset of fields (required for output) of filtered consequences
        results_by_variant[variant_identifier].extend([
            (variant_identifier, c['gene_id'], c['gene_symbol'], most_severe_consequence_term,
             c.get('distance', 0) if report_distance else 0)
            for c in consequences
        ])

    return results_by_variant


def get_variants_without_consequences(results_by_variant):
    """Returns a list of variants for which no consequences were found."""
    return sorted({
        variant_id
        for variant_id, list_of_consequences in results_by_variant.items()
        if len(list_of_consequences) == 0
    })


def process_variants(variants, enable_distant_querying=False, report_distance=False):
    """Given a list of variant IDs, return a list of consequence types (each including Ensembl gene name & ID and a
    functional consequence code) for a given variant.

    Args:
        enable_distant_querying: If set to True, an additional VEP query will be performed for variants for which no
            consequences were found during the first iteration, in an attempt to find distant variant consequences.
        report_distance: Whether to report distance to the nearest gene for upstream and downstream gene variants.
            See extract_consequences() for details.
    """

    # First, we query VEP with default parameters, looking for variants affecting protein coding and miRNA transcripts
    # up to a standard distance (5000 nucleotides either way, which is default for VEP) from the variant.
    results_by_variant = {}
    vep_results = query_vep(variants=variants, search_distance=VEP_SHORT_QUERY_DISTANCE)
    results_by_variant = extract_consequences(vep_results=vep_results, acceptable_biotypes={'protein_coding', 'miRNA'},
                                              only_closest=False, results_by_variant=results_by_variant,
                                              report_distance=report_distance)

    # See if there are variants with no consequences up to the default distance
    variants_without_consequences = get_variants_without_consequences(results_by_variant)
    if variants_without_consequences:
        logger.info('Found {} variant(s) without standard consequences: {}'.format(
            len(variants_without_consequences), '|'.join(variants_without_consequences)))

        if enable_distant_querying:
            logger.info('Attempting to find distant consequences for the remaining variants')

            # If there are, we will now do a second round of querying, this time looking only at protein coding biotypes
            # (vs. miRNA *and* protein coding during the first round) up to a distance of 500,000 bases each way.
            if variants_without_consequences:
                distant_vep_results = query_vep(variants=variants_without_consequences,
                                                search_distance=VEP_LONG_QUERY_DISTANCE)
                extract_consequences(vep_results=distant_vep_results, acceptable_biotypes={'protein_coding'},
                                     only_closest=True, results_by_variant=results_by_variant,
                                     report_distance=report_distance)

            # See if there are still variants with no consequences, even up to a wide search window
            variants_without_consequences = get_variants_without_consequences(results_by_variant)
            if variants_without_consequences:
                logger.info('After distant querying, still remaining {} variant(s) without consequences: {}'.format(
                    len(variants_without_consequences), '|'.join(variants_without_consequences)
                ))

    # Yield all consequences for all variants. Note they are not grouped by variant, all consequences are yielded in a
    # common sequence.
    for variant_id, variant_consequences in results_by_variant.items():
        for consequence_to_yield in deduplicate_list(variant_consequences):
            yield consequence_to_yield


def main():
    # Parse command line arguments
    args = parser.parse_args()

    # Load variants to query from STDIN
    variants_to_query = [colon_based_id_to_vep_id(v) for v in sys.stdin.read().splitlines()]

    # Query VEP with all variants at once (for the purpose of efficiency), print out the consequences to STDOUT.
    consequences = process_variants(variants_to_query,
                                    enable_distant_querying=args.enable_distant_querying,
                                    report_distance=args.report_distance)
    for variant_id, gene_id, gene_symbol, consequence_term, distance in consequences:
        # The second column, set statically to 1, is not used, and is maintained for compatibility purposes
        print('\t'.join([vep_id_to_colon_id(variant_id), '1', gene_id, gene_symbol, consequence_term, str(distance)]))

    logger.info('Successfully processed {} variants'.format(len(variants_to_query)))


if __name__ == '__main__':
    main()
