#!/usr/bin/env python3

import argparse
from collections import Counter
import re
import sys

import eva_cttv_pipeline.clinvar_xml_utils as clinvar_xml_utils

SIG_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, conflicting interpretations': 1,
    'criteria provided, single submitter': 1,
}

parser = argparse.ArgumentParser()
parser.add_argument('--clinvar-xml', required=True)
args = parser.parse_args()


def add_transitions(transitions_counter, transition_chain):
    """Increments the count of a particular flow in Sankey diagram."""
    for transition_from, transition_to in zip(transition_chain, transition_chain[1:]):
        transitions_counter[(transition_from, transition_to)] += 1


def find_attribute(rcv, xpath, attribute_name):
    """Find an attribute in the RCV record which can have either zero or one occurrence. Return a textual representation
    of the attribute, including special representations for the case of zero or multiple, constructed using the
    attribute_name parameter."""

    attributes = rcv.findall(xpath)
    if len(attributes) == 0:
        return '{} missing'.format(attribute_name)
    elif len(attributes) == 1:
        return attributes[0].text
    else:
        return '{} multiple'.format(attribute_name)


def review_status_stars(review_status):
    black_stars = SIG_STARS.get(review_status, 0)
    white_stars = 4 - black_stars
    return '★' * black_stars + '☆' * white_stars


# The dicts store transition counts for the Sankey diagrams. Keys are (from, to), values are transition counts.
# Sankey diagrams can be visualised with SankeyMatic (see http://www.sankeymatic.com/build/).
variant_type_transitions, clin_sig_transitions, review_status_transitions, inheritance_mode_transitions,\
    allele_origin_transitions = Counter(), Counter(), Counter(), Counter(), Counter()
all_clinical_significance_levels = set()


# ClinVar XML have the following top-level structure:
#   <ReleaseSet>
#     <ClinVarSet>...</ClinVarSet>
#     <ClinVarSet>...</ClinVarSet>
#     ...
#   </ReleaseSet>
# To not load everything into memory, we use iterparse, wait for a complete ClinVarSet element, and then remove it from
# the tree because we don't need it anymore.

elements_processed = 0
for rcv in clinvar_xml_utils.iterate_rcv_from_xml(args.clinvar_xml):
    rcv_id = 'RCV{:09}'.format(int(rcv.attrib['ID']))  # FIXME FIXME FIXME WRONG

    # RCV can contain either a MeasureSet, or a GenotypeSet. It must not contain both.
    measure_sets = rcv.findall('MeasureSet')
    genotype_sets = rcv.findall('GenotypeSet')
    if len(measure_sets) == 1 and len(genotype_sets) == 0:
        # Most common case. RCV directly contains one measure set.
        measure_set = measure_sets[0]
        measure_set_type = measure_set.attrib['Type']
        add_transitions(variant_type_transitions, ('RCV', 'MeasureSet', measure_set_type))

        if measure_set_type == 'Variant':
            # Most common case, accounting for >99.95% of all ClinVar records.. Here, we go into details on various
            # attribute distributions.

            # Variant type
            measures = measure_set.findall('Measure')
            assert len(measures) == 1, 'MeasureSet of type Variant must contain exactly one Measure'
            add_transitions(variant_type_transitions, (measure_set_type, measures[0].attrib['Type']))

            # Clinical significance
            clinical_significance = find_attribute(
                rcv, 'ClinicalSignificance/Description', 'ClinicalSignificance')
            all_clinical_significance_levels.add(clinical_significance)
            significance_type = 'Complex' if re.search('[,/]', clinical_significance) else 'Simple'
            add_transitions(clin_sig_transitions, (
                'Variant',
                significance_type,
                clinical_significance,
            ))

            # Review status
            review_status = find_attribute(
                rcv, 'ClinicalSignificance/ReviewStatus', 'ReviewStatus')
            add_transitions(review_status_transitions, (
                'Variant',
                review_status_stars(review_status),
                review_status,
            ))

            # Mode of inheritance
            mode_of_inheritance_xpath = 'AttributeSet/Attribute[@Type="ModeOfInheritance"]'
            mode_of_inheritance = find_attribute(rcv, mode_of_inheritance_xpath, 'ModeOfInheritance')
            if mode_of_inheritance.endswith('multiple'):
                # Having multiple ModeOfInheritance is rare. Log them for further investigation
                all_modes = '|'.join(sorted(mode.text for mode in rcv.findall(mode_of_inheritance_xpath)))
                print(f'Multiple ModeOfInheritance\t{rcv_id}\t{all_modes}')
            add_transitions(inheritance_mode_transitions, (
                'Variant',
                mode_of_inheritance if mode_of_inheritance.endswith('missing') else 'ModeOfInheritance present',
            ))
            if not mode_of_inheritance.endswith('missing'):
                add_transitions(inheritance_mode_transitions, (
                    'ModeOfInheritance present', mode_of_inheritance
                ))

    elif len(measure_sets) == 0 and len(genotype_sets) == 1:
        # RCV directly contains one genotype set.
        genotype_set = genotype_sets[0]
        add_transitions(variant_type_transitions, ('RCV', 'GenotypeSet', genotype_set.attrib['Type']))
    else:
        raise AssertionError('RCV must contain either exactly one measure set, or exactly one genotype set')

    allele_origins = {origin.text for origin in rcv.findall('ObservedIn/Sample/Origin')}
    if len(allele_origins) == 0:
        add_transitions(allele_origin_transitions, ('RCV', 'No allele origin'))
    else:
        allele_origins_count = 'Single allele origin' if len(allele_origins) == 1 else 'Multiple allele origins'
        allele_origins_text = ','.join(sorted(allele_origins))
        add_transitions(allele_origin_transitions, ('RCV', allele_origins_count, allele_origins_text))

    # Track the number of already processed elements
    elements_processed += 1
    if elements_processed % 10000 == 0:
        print('Processed {} elements'.format(elements_processed), file=sys.stderr)

# Output the code for Sankey diagram. Transitions are sorted in decreasing number of counts, so that the most frequent
# cases are on top.
for transitions_counter in (variant_type_transitions, clin_sig_transitions, review_status_transitions,
                            inheritance_mode_transitions, allele_origin_transitions):
    print()
    for (transition_from, transition_to), count in sorted(transitions_counter.items(), key=lambda x: -x[1]):
        print('{transition_from} [{count}] {transition_to}'.format(**locals()))

print('\n\nAll clinical significance levels:')
for clin_sig in sorted(all_clinical_significance_levels):
    print(clin_sig)
