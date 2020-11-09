"""Contains utilities and classes to parse the ClinVar XML and convert the records into internal representation via
ClinVarDataset, ClinVarRecord, and ClinVarRecordMeasure classes."""

import logging
import re
import xml.etree.ElementTree as ElementTree

from eva_cttv_pipeline.file_utils import open_file

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def iterate_rcv_from_xml(clinvar_xml):
    """Iterates through ClinVar XML (possibly gzipped) and yields complete <ReferenceClinVarAssertion> records."""
    with open_file(clinvar_xml, 'rt') as fh:
        for event, elem in ElementTree.iterparse(fh):
            # Wait until we have built a complete ClinVarSet element
            if elem.tag != 'ClinVarSet':
                continue

            # Go to a ReferenceClinVarAssertion element. This corresponds to a single RCV record, the main unit of
            # ClinVar. There should only be one such record per ClinVarSet.
            rcv = find_mandatory_unique_element(elem, 'ReferenceClinVarAssertion')

            # Return the complete record and then remove the processed element from the tree to save memory
            yield rcv
            elem.clear()


def find_elements(node, xpath, allow_zero=True, allow_multiple=True):
    """Attempt to find child elements in a node by xpath. Raise exceptions if conditions are violated. Return a
    (possibly empty) list of elements."""
    all_elements = node.findall(xpath)
    if (len(all_elements) == 0 and not allow_zero) or (len(all_elements) > 1 and not allow_multiple):
        raise AssertionError(f'Found {len(all_elements)} instances of {xpath} in {node}, which is not allowed')
    return all_elements


def find_mandatory_unique_element(node, xpath):
    """Attempt to find a child element by xpath which must have exactly one occurrence, otherwise throw an exception."""
    return find_elements(node, xpath, allow_zero=False, allow_multiple=False)[0]


def find_optional_unique_element(node, xpath):
    """Attempt to find a child element by xpath which must have either 0 or 1 occurrence. Return the element, or None if
    not found."""
    elements = find_elements(node, xpath, allow_zero=True, allow_multiple=False)
    return elements[0] if elements else None


def pubmed_refs_to_urls(pubmed_refs):
    """Convert a list of PubMed identifiers to complete URLs."""
    return [f'http://europepmc.org/abstract/MED/{r}' for r in pubmed_refs]


class ClinVarDataset:
    """Iterate through records in ClinVar XML dump and convert them into internal ClinVarRecord representation."""
    def __init__(self, clinvar_xml):
        self.clinvar_xml = clinvar_xml

    def __iter__(self):
        for rcv in iterate_rcv_from_xml(self.clinvar_xml):
            yield ClinVarRecord(rcv)


class ClinVarRecord:
    """Instances of this class hold data on individual ClinVar records. See also:
    * /clinvar-variant-types/README.md for the in-depth explanation of ClinVar data model;
    * Issue https://github.com/EBIvariation/eva-opentargets/issues/127 for the most recent discussions on changing
      support of different ClinVar record types."""

    # A score for the review status of the assigned clinical significance ranges from 0 to 4 and corresponds to the
    # number of "gold stars" displayed on ClinVar website. See details here:
    # https://www.ncbi.nlm.nih.gov/clinvar/docs/details/#review_status
    score_map = {
        "no assertion provided": 0,
        'no assertion criteria provided': 0,
        'criteria provided, single submitter': 1,
        'criteria provided, conflicting interpretations': 1,
        'criteria provided, multiple submitters, no conflicts': 2,
        'reviewed by expert panel': 3,
        'practice guideline': 4,
    }

    def __init__(self, rcv):
        """Initialise a ClinVar record object from an RCV XML record."""
        self.rcv = rcv

        # Add a list of traits
        self.trait_set = []
        for trait in find_elements(self.rcv, './TraitSet/Trait[@Type="Disease"]'):
            self.trait_set.append(ClinVarTrait(trait, self))

        # We are currently only processing MeasureSets of type Variant which are included directly in the RCV record.
        # Some other options (currently not supported) are:
        # * MeasureSet of types "Haplotype", "Phase unknown", or "Distinct chromosomes"
        # * GenotypeSet, which contains an assertion about a group of variants from different chromosome copies, with
        #   the type of be either a "CompoundHeterozygote" or a "Diplotype"
        variant_measure = find_optional_unique_element(self.rcv, './MeasureSet[@Type="Variant"]/Measure')
        if not variant_measure:
            self.measure = None
        else:
            self.measure = ClinVarRecordMeasure(variant_measure, self)

    def __str__(self):
        return f'ClinVarRecord object with accession {self.accession}'

    @property
    def accession(self):
        return find_mandatory_unique_element(self.rcv, './ClinVarAccession').attrib['Acc']

    @property
    def date(self):
        """This tracks the latest update date, counting even minor technical updates."""
        return self.rcv.attrib['DateLastUpdated']

    @property
    def last_evaluated_date(self):
        """This tracks the latest (re)evaluation date for the clinical interpretation.
        See https://github.com/opentargets/platform/issues/1161#issuecomment-683938510 for details."""
        # The DateLastEvaluated attribute is not always present. In this case, this property will be None.
        return find_mandatory_unique_element(self.rcv, './ClinicalSignificance').attrib.get('DateLastEvaluated')

    @property
    def review_status(self):
        """Return a review status text for the assigned clinical significance. See score_map above for the list of
        possible values."""
        review_status = find_mandatory_unique_element(self.rcv, './ClinicalSignificance/ReviewStatus').text
        assert review_status in self.score_map, f'Unknown review status {review_status} in RCV {self.accession}'
        return review_status

    @property
    def score(self):
        """Return a score (star rating) for the assigned clinical significance. See score_map above."""
        return self.score_map[self.review_status]

    @property
    def mode_of_inheritance(self):
        """Return a (possibly empty) list of modes of inheritance for a given ClinVar record."""
        return sorted({
            elem.text for elem in find_elements(self.rcv, './AttributeSet/Attribute[@Type="ModeOfInheritance"]')
        })

    @property
    def traits(self):
        """Returns a list of traits associated with the ClinVar record, in the form of Trait objects."""
        return self.trait_set

    @property
    def observed_pubmed_refs(self):
        return [int(elem.text)
                for elem in find_elements(self.rcv, './ObservedIn/ObservedData/Citation/ID[@Source="PubMed"]')]

    @property
    def clinical_significance_raw(self):
        """The original clinical significance string as stored in ClinVar. Example: 'Benign/Likely benign'."""
        return find_mandatory_unique_element(self.rcv, './ClinicalSignificance/Description').text

    @property
    def clinical_significance_list(self):
        """The normalised list of all clinical significance values. The original value is (1) split into multiple values
        by two delimiters: ('/', ', '), (2) converted into lowercase and (3) sorted lexicographically. Example:
        'Benign/Likely benign, risk_factor' → ['benign', 'likely benign', 'risk factor']. See /clinvar-variant-types/
        README.md for further explanation."""
        return sorted(re.split('/|, ', self.clinical_significance_raw.lower().replace('_', ' ')))

    @property
    def allele_origins(self):
        return {elem.text for elem in find_elements(self.rcv, './ObservedIn/Sample/Origin')}


class ClinVarTrait:
    """Represents a single ClinVar trait (usually a disease), with the corresponding database and Pubmed
    cross-references."""

    def __init__(self, trait_xml, clinvar_record):
        self.trait_xml = trait_xml
        self.clinvar_record = clinvar_record

    def __str__(self):
        return f'ClinVarTrait object with name {self.name} from ClinVar record {self.clinvar_record.accession}'

    @property
    def name(self):
        return find_mandatory_unique_element(self.trait_xml, './Name/ElementValue[@Type="Preferred"]').text

    @property
    def pubmed_refs(self):
        """Trait-specific PubMed references, contained inside a Trait entity."""
        return [int(elem.text) for elem in find_elements(self.trait_xml, './Citation/ID[@Source="PubMed"]')]

    @property
    def xrefs(self):
        return [(elem.attrib['DB'], elem.attrib['ID'], elem.attrib.get('Status', 'current').lower())
                for elem in find_elements(self.trait_xml, './XRef')]


class ClinVarRecordMeasure:
    """This class represents individual ClinVar record "measures". Measures are essentially isolated variants, which can
    be combined into either MeasureSets (include one or more Measures) or GenotypeSets. For a detailed description of
    ClinVar data model, see /clinvar-variant-types/."""

    def __init__(self, measure_xml, clinvar_record):
        self.measure_xml = measure_xml
        self.clinvar_record = clinvar_record

    @property
    def name(self):
        return find_mandatory_unique_element(self.measure_xml, './Name/ElementValue[@Type="Preferred"]').text

    @property
    def preferred_gene_symbols(self):
        return [elem.text for elem in find_elements(
            self.measure_xml, './MeasureRelationship/Symbol/ElementValue[@Type="Preferred"]')]

    @property
    def hgnc_ids(self):
        return [elem.attrib['ID'] for elem in find_elements(
            self.measure_xml, './MeasureRelationship/XRef[@DB="HGNC"]')]

    @property
    def rs_id(self):
        """Returns dbSNP rsID found for the record measure."""
        rs_ids = ['rs' + elem.attrib['ID'] for elem in find_elements(self.measure_xml, './XRef[@DB="dbSNP"]')]
        if len(rs_ids) == 0:
            return None
        elif len(rs_ids) == 1:
            return rs_ids[0]
        else:
            logger.warning(f'Found multiple RS IDs for {self.clinvar_record}, this is not yet supported')
            return None

    @property
    def nsv_id(self):
        nsv_ids = [elem.attrib['ID']
                   for elem in find_elements(self.measure_xml, './XRef[@DB="dbVar"]')
                   if elem.attrib['ID'].startswith('nsv')]
        if len(nsv_ids) == 0:
            return None
        elif len(nsv_ids) == 1:
            return nsv_ids[0]
        else:
            logger.warning(f'Found multiple NSV IDs for {self.clinvar_record}, this is not yet supported')
            return None

    @property
    def hgvs(self):
        return [
            elem.text
            for elem in find_elements(self.measure_xml, './AttributeSet/Attribute')
            if elem.attrib['Type'].startswith('HGVS')
        ]

    @property
    def variant_type(self):
        return self.measure_xml.attrib['Type']

    @property
    def pubmed_refs(self):
        """Variant-specific PubMed references, contained inside a Measure entity."""
        return [int(elem.text) for elem in find_elements(self.measure_xml, './Citation/ID[@Source="PubMed"]')]

    @property
    def chr(self):
        return self.sequence_location_helper('Chr')

    @property
    def vcf_pos(self):
        return self.sequence_location_helper('positionVCF')

    @property
    def vcf_ref(self):
        return self.sequence_location_helper('referenceAlleleVCF')

    @property
    def vcf_alt(self):
        return self.sequence_location_helper('alternateAlleleVCF')

    @property
    def has_complete_coordinates(self):
        return self.chr and self.vcf_pos and self.vcf_ref and self.vcf_alt

    def sequence_location_helper(self, attr):
        if self.variant_type == 'Translocation':
            # TODO: Translocations have multiple locations and are not supported.
            # TODO: https://github.com/EBIvariation/eva-opentargets/issues/171
            return None
        sequence_locations = find_elements(self.measure_xml, './SequenceLocation[@Assembly="GRCh38"]')
        if len(sequence_locations) != 1:
            # TODO: Support variants with multiple locations (for example, chrX/chrY).
            # TODO: https://github.com/EBIvariation/eva-opentargets/issues/172
            return None
        return sequence_locations[0].attrib.get(attr)
