"""
The regular expressions in this module are used to parse ClinVar's identifiers for describing repeat expansion variants.
The pyhgvs module cannot be applied here because not all expression used by ClinVar are actually valid HGVS, and that
module imposes strict validation. Hence, custom regular expressions are necessary.
"""

import logging
import re

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# Common part for all HGVS-like transcript definitions, e.g. 'NM_001256054.2(C9orf72):'
hgvs_like_transcript_part = (
    r'(?P<transcript_id>[a-zA-Z][a-zA-Z0-9_]+)'  # Transcript accession                      NM_001256054
    r'\.'                                        # Delimiter, transcript accession/version   .
    r'[0-9]+'                                    # Transcript version                        2

    r'(?:'                                       # Non-capturing group for the gene symbol
    r'\('                                            # Opening parenthesis                   (
    r'[a-zA-Z0-9_.]+'                                # Gene name                             C9orf72
    r'\)'                                            # Closing parenthesis                   )
    r')?'                                        # The entire gene symbol part is optional
    
    r':'                                         # Delimiter, transcript/variant info        :
)

# Common part for start and end coordinate pivots. Pivots for introns and other cases where a coordinate in a noncoding
# region of mRNA needs to be addressed relative to the coding regions. For example, c.87 means the 87th coding base of
# mRNA (not counting introns). In comparison, c.87+14 means that base number 87 is the last base of a particular exon,
# and the base addressed by the coordinate is 14 bases downstream of the pivot base. In this case, an intron repeat
# expansion variant might be addressed as c.87+14_c.87+17. In this case, we don't want the pivots (87), but only the
# actual coordinates (+14 and +17). To do that, we have a regular expression which captures the pivot part.
coordinate_pivot_part = (
    r'(?:'       # Non-capturing group for coordinate pivot part
    r'[-+]?'     # Coordinate pivot may start with a minus or plus
    r'[0-9]+'    # Then it contains a number
    r'(?=[-+])'  # It then must be followed with either a plus or minus.  This is a positive lookahead and is *not*
                 # part of the coordinate pivot, hence the (?=...) notation
    r')?'        # The coordinate pivot is always optional
)

# Pattern 1. HGVS-like notation for coding or genomic coordinates. This pattern is used for most variants. From it we
# can always extract transcript ID and start coord, and sometimes also end coord and repeat unit sequence.
# Example: 'NM_001256054.2(C9orf72):c.-45+63_-45+80GGGGCC(2_25)'
re_hgvs_like_transcript_or_genomic = re.compile(
    hgvs_like_transcript_part +      # Transcript definition                                NM_001256054.2(C9orf72):
    r'[gc]'                          # Coordinate type, genomic or coding                   c
    r'\.'                            # Delimiter, coordinate type/coordinate                .

    + coordinate_pivot_part +        # Start coordinate pivot, optional                     -45
    r'\*?'                           # Sometimes there is an asterisk in front of the
                                     # coordinate (special case, not important)
    r'(?P<start_coord>[+-]?[0-9]+)'  # Start coordinate                                     +63

    r'(?:'                           # Non-capruting group for end coordinate
    r'_'                                 # Delimiter: start coordinate/end coordinate       _
    + coordinate_pivot_part +            # End coordinate pivot, optional                   -45
    r'\*?'                               # Coordinate may start with an asterisk
    r'(?P<end_coord>[+-]?[0-9]+)'        # End coordinate                                   +80
    r')?'                            # The entire end coordinate part is optional

    r'(?P<sequence>[{}]*)'.format(   # Repeat unit sequence, optional                       GGGGCC or
        IUPACAmbiguousDNA.letters)   # IUPAC ambiguity codes are supported                  CGN
)

# Pattern 2. HGVS-like notation for protein coordinates. This pattern is used for a few variants. From it we do not
# extract any useful information, but always assume that such a notation signified a trinucleotide repeat type.
# Example: NP_002964.3:p.Gln166(>=33)
re_hgvs_like_protein = re.compile(hgvs_like_transcript_part + r'p\.')

# Pattern 3. Human-readable description notation. This pattern is used for about a dozen variants. From it we only
# extract the repeat unit sequence. Example: 'ATXN8, (CAG)n REPEAT EXPANSION'
re_description = re.compile(
    r'\('
    r'(?P<sequence>[{}]+)'.format(IUPACAmbiguousDNA.letters) +
    r'\)n'
    r'(?: REPEAT)? EXPANSION'
)


def parse_variant_identifier(variant_name):
    """
    Parses variant identifier and extract certain characteristics:
        * TranscriptID: NCBI RefSeq transcript ID, e.g. NM_000044
        * CoordinateSpan: the distance between start and end coordinates described in the HGVS-like notation.
              E. g., for 'NM_000044.4(AR):c.172_174CAG(7_34) (p.Gln66_Gln80del)', it will be 174 - 172 + 1 = 3
        * RepeatUnitLength: length of the sequence being repeated.
              E. g., for 'NC_000004.11:g.3076606GCA[27_35]', it will be len('GCA') = 3
        * IsProteinHGVS: this field simply reflects whether the variant was defined using a protein HGVS-like notation.
              E.g., for 'NP_002964.3:p.Gln166(>=33)' it will be TRUE, and for the two examples above FALSE.
    """
    # A variant "name", or identifier, can come from either of three patterns described in the section above.
    transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs = None, None, None, False

    # Try to match HGVS-like transcript/genomic ID
    match = re_hgvs_like_transcript_or_genomic.search(variant_name)
    if match:
        transcript_id, start_coord, end_coord, sequence = \
            match.group('transcript_id'), match.group('start_coord'), match.group('end_coord'), match.group('sequence')
        if transcript_id.startswith('NM'):  # We are only interested in RefSeq mRNA transcripts for querying
            transcript_id = match.group('transcript_id')
        if start_coord and end_coord:  # If both start *and* end coordinates are present, we can calculate the span
            coordinate_span = int(end_coord) - int(start_coord) + 1
        if sequence:  # If the repeat unit is present, we can calculate its length directly
            repeat_unit_length = len(sequence)
        return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs

    # Is this a protein HGVS-like notation?
    if re_hgvs_like_protein.search(variant_name):
        is_protein_hgvs = True
        return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs

    # Is this a human-readable description?
    match = re_description.search(variant_name)
    if match:
        repeat_unit_length = len(match.group('sequence'))
        return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs

    logger.warning('ClinVar identifier did not match any of the regular expressions: {}'.format(variant_name))
    return transcript_id, coordinate_span, repeat_unit_length, is_protein_hgvs
