from collections import Counter

from eva_cttv_pipeline import clinvar_xml_utils
from eva_cttv_pipeline.trait_mapping.trait import Trait


def parse_trait_names(filepath: str) -> list:
    """For a file containing ClinVar records in the XML format, return a list of Traits for the records in the file.
    Each Trait object contains trait name, how many times it occurs in the input file, and whether it is linked to an NT
    expansion variant.

    Trait occurrence count is calculated based on all unique (RCV, trait name) tuples in the input file. This is because
    each such tuple will, generally speaking, correspond to one output evidence string. So if we want to gauge which
    trait names are more important to curate, we need to consider how many such tuples it appears in.

    Traits which are implicated in "Microsatellite" variants are marked using a special field, because a subset of
    microsatellites are NT expansion variants, and their curation is of highest importance even if the number of records
    which they are linked to is low.

    :param filepath: Path to a gzipped file containing ClinVar XML dump.
    :return: A list of Trait objects."""

    # Tracks how many times a trait name occurs in ClinVar
    trait_name_counter = Counter()

    # Tracks all traits which are at least once implicated in "NT expansion", or nucleotide repeat expansion, variants.
    # Their curation is of highest importance regardless of how many records they are actually associated with.
    nt_expansion_traits = set()

    for clinvar_record in clinvar_xml_utils.ClinVarDataset(filepath):
        trait_names = set(trait.preferred_or_other_valid_name.lower()
                          for trait in clinvar_record.traits_with_valid_names)
        for trait_name in trait_names:
            trait_name_counter[trait_name] += 1
        if clinvar_record.measure and clinvar_record.measure.is_repeat_expansion_variant:
            nt_expansion_traits |= trait_names

    # Count trait occurrences
    traits = []
    for trait_name, trait_frequency in trait_name_counter.items():
        if trait_name == '-':
            print('Skipped {} missing trait names'.format(trait_frequency))
            continue
        associated_with_nt_expansion = trait_name in nt_expansion_traits
        traits.append(Trait(name=trait_name, frequency=trait_frequency,
                            associated_with_nt_expansion=associated_with_nt_expansion))

    return traits
