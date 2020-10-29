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

    # Tracks unique (RCV, trait name) tuples
    unique_association_tuples = set()

    # Tracks all traits which are at least once implicated in "NT expansion", or nucleotide repeat expansion, variants.
    # Their curation is of highest importance regardless of how many records they are actually associated with.
    nt_expansion_traits = set()

    for clinvar_record in clinvar_xml_utils.ClinVarDataset(filepath):
        traits = set(trait.name for trait in clinvar_record.traits)
        unique_association_tuples |= {(clinvar_record.accession, trait) for trait in traits}
        # FIXME not all microsatellites are actually NT expansions!
        if clinvar_record.measure is not None and clinvar_record.measure.variant_type == 'Microsatellite':
            nt_expansion_traits |= traits

    # Count trait occurrences
    trait_names = [t[1] for t in unique_association_tuples]
    traits = []
    for trait_name, trait_frequency in Counter(trait_names).items():
        if trait_name == '-':
            print('Skipped {} missing trait names'.format(trait_frequency))
            continue
        associated_with_nt_expansion = trait_name in nt_expansion_traits
        traits.append(Trait(name=trait_name.lower(), frequency=trait_frequency,
                            associated_with_nt_expansion=associated_with_nt_expansion))

    return traits
