from datetime import datetime
import unittest
from types import SimpleNamespace

from eva_cttv_pipeline import clinvar, efo_term, clinvar_to_evidence_strings, \
    evidence_strings, consequence_type

import tests.test_config as test_config
from tests import test_clinvar_to_evidence_strings

DATE_LAST_UPDATED = 1412982000000
DATE_ASSERTED = datetime.fromtimestamp((DATE_LAST_UPDATED / 1000)).isoformat()


MAPPINGS = SimpleNamespace()


def _get_test_cellbase_record_gene():
    return {'end': 23866396, 'clinvarSet': {'referenceClinVarAssertion': {'dateLastUpdated': DATE_LAST_UPDATED, 'dateCreated': 1405551600000, 'id': 300575, 'traitSet': {'trait': [{'xref': [{'id': 'C0018798', 'db': 'MedGen', 'status': 'CURRENT'}, {'id': '140500', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}, {'id': '234750', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'name': [{'elementValue': {'value': 'Malformation of the heart', 'type': 'Preferred'}}, {'xref': [{'id': '140500', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}, {'id': '234750', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'elementValue': {'value': 'HEART, MALFORMATION OF', 'type': 'Alternate'}}, {'elementValue': {'value': 'Congenital heart defect', 'type': 'Alternate'}}], 'id': 15882, 'type': 'Disease'}], 'id': 16428, 'type': 'Disease'}, 'measureSet': {'name': [{'elementValue': {'value': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser)', 'type': 'preferred name'}}], 'id': 139663, 'measure': [{'id': 150127, 'cytogeneticLocation': ['14q11.2'], 'measureRelationship': [{'symbol': [{'elementValue': {'value': 'MYH6', 'type': 'Preferred'}}], 'xref': [{'id': '4624', 'db': 'Gene', 'status': 'CURRENT'}, {'id': '160710', 'db': 'OMIM', 'status': 'CURRENT', 'type': 'MIM'}], 'type': 'variant in gene', 'name': [{'elementValue': {'value': 'myosin, heavy chain 6, cardiac muscle, alpha', 'type': 'Preferred'}}], 'sequenceLocation': [{'stop': 23877485, 'accession': 'NC_000014.8', 'assembly': 'GRCh37', 'chr': '14', 'strand': '-', 'start': 23851198}, {'stop': 23409621, 'accession': 'NC_000014.9', 'assembly': 'GRCh38', 'chr': '14', 'strand': '-', 'start': 23380732}]}], 'sequenceLocation': [{'stop': 23866396, 'accession': 'NC_000014.8', 'variantLength': 1, 'assembly': 'GRCh37', 'referenceAllele': 'T', 'chr': '14', 'start': 23866396, 'alternateAllele': 'C'}, {'stop': 23397187, 'accession': 'NC_000014.9', 'variantLength': 1, 'assembly': 'GRCh38', 'referenceAllele': 'T', 'chr': '14', 'start': 23397187, 'alternateAllele': 'C'}], 'attributeSet': [{'attribute': {'change': 'c.2033A>G', 'value': 'LRG_389t1:c.2033A>G', 'type': 'HGVS, coding, LRG'}}, {'attribute': {'change': 'c.2033A>G', 'value': 'NM_002471.3:c.2033A>G', 'type': 'HGVS, coding, RefSeq'}}, {'attribute': {'change': 'g.16091A>G', 'value': 'LRG_389:g.16091A>G', 'type': 'HGVS, genomic, LRG'}}, {'attribute': {'change': 'g.16091A>G', 'value': 'NG_023444.1:g.16091A>G', 'type': 'HGVS, genomic, RefSeqGene'}}, {'attribute': {'change': 'g.23397187T>C', 'value': 'NC_000014.9:g.23397187T>C', 'integerValue': 38, 'type': 'HGVS, genomic, top level'}}, {'attribute': {'change': 'g.23866396T>C', 'value': 'NC_000014.8:g.23866396T>C', 'integerValue': 37, 'type': 'HGVS, genomic, top level, previous'}}, {'attribute': {'change': 'p.Asn678Ser', 'value': 'LRG_389p1:p.Asn678Ser', 'type': 'HGVS, protein'}}, {'attribute': {'change': 'p.Asn678Ser', 'value': 'NP_002462.2:p.Asn678Ser', 'type': 'HGVS, protein, RefSeq'}}, {'xref': [{'id': 'SO:0001583', 'db': 'Sequence Ontology', 'status': 'CURRENT'}, {'id': 'NM_002471.3:c.2033A>G', 'db': 'RefSeq', 'status': 'CURRENT'}], 'attribute': {'value': 'missense variant', 'type': 'MolecularConsequence'}}, {'attribute': {'value': 'N678S', 'type': 'ProteinChange1LetterCode'}}], 'xref': [{'id': '515726230', 'db': 'dbSNP', 'status': 'CURRENT', 'type': 'rs'}], 'name': [{'elementValue': {'value': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser)', 'type': 'Preferred'}}], 'type': 'single nucleotide variant'}], 'type': 'Variant'}, 'assertion': {'type': 'VARIATION_TO_DISEASE'}, 'clinicalSignificance': {'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER', 'description': 'Likely pathogenic'}, 'observedIn': [{'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'taxonomyId': 9606, 'value': 'human'}, 'origin': 'germline', 'numberTested': 1}, 'observedData': [{'id': 3592924, 'attribute': {'value': 'not provided', 'type': 'Description'}}]}, {'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'taxonomyId': 9606, 'value': 'human'}, 'origin': 'inherited', 'numberTested': 1}, 'observedData': [{'id': 3642470, 'attribute': {'value': 'not provided', 'type': 'Description'}}]}], 'clinVarAccession': {'acc': 'RCV000128628', 'dateUpdated': DATE_LAST_UPDATED, 'version': 1, 'type': 'RCV'}, 'recordStatus': 'current'}, 'title': 'NM_002471.3(MYH6):c.2033A>G (p.Asn678Ser) AND Malformation of the heart', 'recordStatus': 'current', 'id': 3829163, 'clinVarAssertion': [{'id': 299693, 'traitSet': {'trait': [{'xref': [{'id': 'C14.240.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}, {'xref': [{'id': 'C14.280.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}, {'xref': [{'id': 'C16.131.240.400', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'name': [{'elementValue': {'value': 'not provided', 'type': 'Preferred'}}], 'type': 'Disease'}], 'type': 'Disease'}, 'measureSet': {'measure': [{'attributeSet': [{'attribute': {'value': 'NM_002471.3:c.2033A>G', 'type': 'HGVS'}}], 'type': 'Variation'}], 'type': 'Variant'}, 'comment': [{'value': 'this variant was identified in a 4 generation family with multiple members affected with congenital heart defects (multiple types)'}], 'assertion': {'type': 'variation to disease'}, 'submissionName': 'cardiac_targeted_resequencing', 'clinVarSubmissionID': {'submitterDate': 1405465200000, 'submitter': 'Laboratory for Genetics of Human Development Center for Human Genetics, Catholic University of Leuven', 'localKey': 'NM_002471.3:c.2033A>G|LABGENHUDEVKULEUVEN'}, 'observedIn': [{'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'value': 'human'}, 'origin': 'germline', 'numberTested': 1}, 'observedData': [{'attribute': {'value': 'not provided', 'type': 'Description'}}]}, {'method': [{'methodType': 'NOT_PROVIDED'}], 'sample': {'affectedStatus': 'not provided', 'species': {'value': 'human'}, 'origin': 'inherited', 'numberTested': 1}, 'observedData': [{'attribute': {'value': 'not provided', 'type': 'Description'}}]}], 'recordStatus': 'current', 'clinVarAccession': {'orgID': 500109, 'acc': 'SCV000172246', 'dateUpdated': 1405638000000, 'version': 1, 'type': 'SCV'}, 'clinicalSignificance': {'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER', 'comment': [{'value': 'Converted during submission to Likely pathogenic.'}], 'description': ['probable-pathogenic']}}]}, 'reference': 'T', 'chromosome': '14', 'alternate': 'C', 'start': 23866396, 'annot': {'end': 23866396, 'hgvs': ['ENST00000405093.3:c.2033A>G', 'ENSP00000386041.3:p.Asn678Ser', 'ENST00000356287.3:c.2033A>G', 'ENSP00000348634.3:p.Asn678Ser'], 'alternativeAllele': 'C', 'chromosome': '14', 'referenceAllele': 'T', 'start': 23866396, 'consequenceTypes': [{'ensemblTranscriptId': 'ENST00000405093', 'biotype': 'protein_coding', 'soTerms': [{'soName': 'missense_variant', 'soAccession': 'SO:0001583'}], 'geneName': 'MYH6', 'codon': 'aAt/aGt', 'proteinSubstitutionScores': [{'description': 'deleterious', 'source': 'Sift', 'score': 0.0}, {'description': 'probably_damaging', 'source': 'Polyphen', 'score': 0.999}], 'cDnaPosition': 2104, 'aaPosition': 678, 'aaChange': 'N/S', 'cdsPosition': 2033, 'ensemblGeneId': 'ENSG00000197616', 'strand': '-'}, {'ensemblTranscriptId': 'ENST00000356287', 'biotype': 'protein_coding', 'soTerms': [{'soName': 'missense_variant', 'soAccession': 'SO:0001583'}], 'geneName': 'MYH6', 'codon': 'aAt/aGt', 'proteinSubstitutionScores': [{'description': 'deleterious', 'source': 'Sift', 'score': 0.0}, {'description': 'probably_damaging', 'source': 'Polyphen', 'score': 0.999}], 'cDnaPosition': 2063, 'aaPosition': 678, 'aaChange': 'N/S', 'cdsPosition': 2033, 'ensemblGeneId': 'ENSG00000197616', 'strand': '-'}]}}


def get_args_CTTVGeneticsEvidenceString_init():
    cellbase_record = _get_test_cellbase_record_gene()

    clinvarRecord = \
        clinvar.ClinvarRecord(mappings=test_clinvar_to_evidence_strings.MAPPINGS,
                                     a_dictionary=cellbase_record['clinvarSet'])

    report = clinvar_to_evidence_strings.Report()

    trait = SimpleNamespace()
    trait.trait_counter = 0
    trait.clinvar_trait_list = [[]]
    trait.efo_list = ['http://www.orpha.net/ORDO/Orphanet_88991']

    ensembl_gene_id = "ENSG00000197616"

    test_args_1 = [clinvarRecord, report, trait, ensembl_gene_id, cellbase_record]

    return test_args_1


# TODO look into why these failed with travis
class CTTVGeneticsEvidenceStringInitTest(unittest.TestCase):
    maxDiff = None
    def setUp(self):
        test_args = get_args_CTTVGeneticsEvidenceString_init()
        self.evidence_string = evidence_strings.CTTVGeneticsEvidenceString(*test_args)

    def test_evidence_string(self):
        test_dict = {'validated_against_schema_version': '1.2.3', 'variant': {'type': 'snp single', 'id': ['http://identifiers.org/dbsnp/rs515726230']}, 'type': 'genetic_association', 'unique_association_fields': {'clinvarAccession': 'RCV000128628', 'alleleOrigin': 'germline', 'gene': 'ENSG00000197616', 'phenotype': 'http://www.orpha.net/ORDO/Orphanet_88991'}, 'evidence': {'variant2disease': {'provenance_type': {'database': {'dbxref': {'id': 'http://identifiers.org/clinvar', 'url': 'http://identifiers.org/clinvar.record/RCV000128628', 'version': '2015-04'}, 'id': 'EVA', 'version': '1.0'}, 'expert': {'status': True, 'statement': 'Primary submitter of data'}}, 'unique_experiment_reference': 'http://europepmc.org/abstract/MED/0', 'is_associated': True, 'evidence_codes': ['http://purl.obolibrary.org/obo/ECO_0000205'], 'urls': [{'url': 'http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628', 'nice_name': 'Further details in ClinVar database'}], 'resource_score': {'type': 'pvalue', 'value': 1e-07, 'method': {'description': 'Not provided by data supplier', 'url': ''}}, 'date_asserted': DATE_ASSERTED}, 'gene2variant': {'provenance_type': {'database': {'dbxref': {'id': 'http://identifiers.org/clinvar', 'url': 'http://identifiers.org/clinvar.record/RCV000128628', 'version': '2015-04'}, 'id': 'EVA', 'version': '1.0'}, 'expert': {'status': True, 'statement': 'Primary submitter of data'}}, 'is_associated': True, 'evidence_codes': ['http://identifiers.org/eco/cttv_mapping_pipeline'], 'functional_consequence': 'http://purl.obolibrary.org/obo/SO_0001583', 'urls': [{'url': 'http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628', 'nice_name': 'Further details in ClinVar database'}], 'date_asserted': DATE_ASSERTED}}, 'access_level': 'public', 'disease': {'id': ['http://www.orpha.net/ORDO/Orphanet_88991']}, 'target': {'id': ['http://identifiers.org/ensembl/ENSG00000197616'], 'target_type': 'http://identifiers.org/cttv.target/gene_variant', 'activity': 'http://identifiers.org/cttv.activity/predicted_damaging'}, 'sourceID': 'eva'}

        test_ev_string = evidence_strings.CTTVEvidenceString(test_dict)

        self.assertEqual(self.evidence_string, test_ev_string)


def get_args_CTTVSomaticEvidenceString_init():
    clinvarRecord = clinvar.ClinvarRecord(mappings=test_clinvar_to_evidence_strings.MAPPINGS,
                                                 a_dictionary={'title': 'NM_031157.2(HNRNPA1):c.1054C>T (p.Arg352Ter) AND Chronic progressive multiple sclerosis', 'clinVarAssertion': [{'measureSet': {'type': 'Variant', 'measure': [{'attributeSet': [{'attribute': {'type': 'HGVS', 'value': 'NM_002136.2:c.898C>T'}}], 'type': 'Variation'}]}, 'id': 286280, 'traitSet': {'type': 'Disease', 'trait': [{'name': [{'elementValue': {'type': 'Preferred', 'value': 'not provided'}}], 'xref': [{'id': 'C10.114.375.500.200', 'db': 'MESH', 'status': 'UNDER_REVIEW'}], 'type': 'Disease'}]}, 'submissionName': 'HNRNPA1_VARIATION', 'clinVarAccession': {'orgID': 504916, 'type': 'SCV', 'dateUpdated': 1409698800000, 'acc': 'SCV000154960', 'version': 1}, 'assertion': {'type': 'variation to disease'}, 'clinVarSubmissionID': {'submitter': 'Demyelinating Disease Laboratories, VA Medical Center and University of Tennessee', 'submitterDate': 1402268400000, 'localKey': 'NP_002127.1:c.898C>T'}, 'clinicalSignificance': {'description': ['probable-pathogenic'], 'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER', 'comment': [{'value': 'Converted during submission to Likely pathogenic.'}]}, 'recordStatus': 'current', 'observedIn': [{'observedData': [{'attribute': {'type': 'Description', 'value': 'not provided'}}], 'sample': {'affectedStatus': 'not provided', 'species': {'value': 'human'}, 'origin': 'somatic', 'numberTested': 1}, 'method': [{'methodType': 'NOT_PROVIDED'}]}]}], 'referenceClinVarAssertion': {'measureSet': {'name': [{'elementValue': {'type': 'preferred name', 'value': 'NM_031157.2(HNRNPA1):c.1054C>T (p.Arg352Ter)'}}], 'measure': [{'attributeSet': [{'attribute': {'type': 'HGVS, coding, RefSeq', 'change': 'c.1054C>T', 'value': 'NM_031157.2:c.1054C>T'}}, {'attribute': {'type': 'HGVS, coding, RefSeq', 'change': 'c.898C>T', 'value': 'NM_002136.2:c.898C>T'}}, {'attribute': {'type': 'HGVS, genomic, RefSeqGene', 'change': 'g.8255C>T', 'value': 'NG_033830.1:g.8255C>T'}}, {'attribute': {'type': 'HGVS, genomic, top level', 'integerValue': 38, 'change': 'g.54283958C>T', 'value': 'NC_000012.12:g.54283958C>T'}}, {'attribute': {'type': 'HGVS, genomic, top level, previous', 'integerValue': 37, 'change': 'g.54677742C>T', 'value': 'NC_000012.11:g.54677742C>T'}}, {'attribute': {'type': 'HGVS, protein, RefSeq', 'change': 'p.Arg300Ter', 'value': 'NP_002127.1:p.Arg300Ter'}}, {'attribute': {'type': 'HGVS, protein, RefSeq', 'change': 'p.Arg352Ter', 'value': 'NP_112420.1:p.Arg352Ter'}}, {'xref': [{'id': 'SO:0001587', 'db': 'Sequence Ontology', 'status': 'CURRENT'}, {'id': 'NM_031157.2:c.1054C>T', 'db': 'RefSeq', 'status': 'CURRENT'}], 'attribute': {'type': 'MolecularConsequence', 'value': 'nonsense'}}, {'attribute': {'type': 'ProteinChange1LetterCode', 'value': 'R352*'}}, {'attribute': {'type': 'ProteinChange1LetterCode', 'value': 'R300*'}}], 'cytogeneticLocation': ['12q13.13'], 'sequenceLocation': [{'variantLength': 1, 'chr': '12', 'start': 54677742, 'assembly': 'GRCh37', 'stop': 54677742, 'alternateAllele': 'T', 'referenceAllele': 'C', 'accession': 'NC_000012.11'}, {'variantLength': 1, 'chr': '12', 'start': 54283958, 'assembly': 'GRCh38', 'stop': 54283958, 'alternateAllele': 'T', 'referenceAllele': 'C', 'accession': 'NC_000012.12'}], 'name': [{'elementValue': {'type': 'Preferred', 'value': 'NM_031157.2(HNRNPA1):c.1054C>T (p.Arg352Ter)'}}], 'xref': [{'type': 'rs', 'id': '483353037', 'db': 'dbSNP', 'status': 'CURRENT'}], 'type': 'single nucleotide variant', 'measureRelationship': [{'symbol': [{'elementValue': {'type': 'Preferred', 'value': 'HNRNPA1'}}], 'name': [{'elementValue': {'type': 'Preferred', 'value': 'heterogeneous nuclear ribonucleoprotein A1'}}], 'xref': [{'id': '3178', 'db': 'Gene', 'status': 'CURRENT'}, {'type': 'MIM', 'id': '164017', 'db': 'OMIM', 'status': 'CURRENT'}], 'type': 'variant in gene', 'sequenceLocation': [{'chr': '12', 'start': 54674487, 'assembly': 'GRCh37', 'stop': 54679029, 'strand': '+', 'accession': 'NC_000012.11'}, {'chr': '12', 'start': 54280695, 'assembly': 'GRCh38', 'stop': 54287086, 'strand': '+', 'accession': 'NC_000012.12'}]}], 'id': 139332}], 'type': 'Variant', 'id': 135606}, 'dateCreated': 1402354800000, 'traitSet': {'type': 'Disease', 'trait': [{'name': [{'elementValue': {'type': 'Preferred', 'value': 'Chronic progressive multiple sclerosis'}, 'xref': [{'id': '230373008', 'db': 'SNOMED CT', 'status': 'CURRENT'}]}], 'xref': [{'id': 'C0393665', 'db': 'MedGen', 'status': 'CURRENT'}], 'type': 'Disease', 'id': 18795}], 'id': 13814}, 'dateLastUpdated': DATE_LAST_UPDATED, 'clinVarAccession': {'type': 'RCV', 'dateUpdated': DATE_LAST_UPDATED, 'acc': 'RCV000122455', 'version': 1}, 'assertion': {'type': 'VARIATION_TO_DISEASE'}, 'clinicalSignificance': {'description': 'Likely pathogenic', 'reviewStatus': 'CLASSIFIED_BY_SINGLE_SUBMITTER'}, 'recordStatus': 'current', 'observedIn': [{'observedData': [{'id': 3574852, 'attribute': {'type': 'Description', 'value': 'not provided'}}], 'sample': {'affectedStatus': 'not provided', 'species': {'taxonomyId': 9606, 'value': 'human'}, 'origin': 'somatic', 'numberTested': 1}, 'method': [{'methodType': 'NOT_PROVIDED'}]}], 'id': 286319}, 'id': 3823230, 'recordStatus': 'current'})

    report = clinvar_to_evidence_strings.Report()

    trait = SimpleNamespace()
    trait.trait_counter = 0
    trait.clinvar_trait_list = [[]]
    trait.efo_list = ['http://www.ebi.ac.uk/efo/EFO_0003840']

    ensembl_gene_id = "ENSG00000135486"

    test_args_1 = (clinvarRecord, report, trait, ensembl_gene_id)

    return test_args_1


class CTTVSomaticEvidenceStringInitTest(unittest.TestCase):
    def setUp(self):
        test_args = get_args_CTTVSomaticEvidenceString_init()
        self.evidence_string = evidence_strings.CTTVSomaticEvidenceString(*test_args)

    def test_evidence_string(self):
        test_dict = {'evidence': {'provenance_type': {'expert': {'status': True, 'statement': 'Primary submitter of data'}, 'database': {'dbxref': {'id': 'http://identifiers.org/clinvar', 'url': 'http://identifiers.org/clinvar.record/RCV000122455', 'version': '2015-04'}, 'id': 'EVA', 'version': '1.0'}}, 'resource_score': {'value': 1, 'type': 'probability'}, 'date_asserted': DATE_ASSERTED, 'is_associated': True, 'known_mutations': [{'functional_consequence': 'http://purl.obolibrary.org/obo/SO_0001587', 'preferred_name': 'stop_gained'}], 'urls': [{'url': 'http://www.ncbi.nlm.nih.gov/clinvar/RCV000122455', 'nice_name': 'Further details in ClinVar database'}], 'evidence_codes': ['http://purl.obolibrary.org/obo/ECO_0000205']}, 'type': 'somatic_mutation', 'sourceID': 'eva_somatic', 'access_level': 'public', 'disease': {'id': ['http://www.ebi.ac.uk/efo/EFO_0003840']}, 'validated_against_schema_version': '1.2.3', 'target': {'activity': 'http://identifiers.org/cttv.activity/predicted_damaging', 'id': ['http://identifiers.org/ensembl/ENSG00000135486'], 'target_type': 'http://identifiers.org/cttv.target/gene_variant'}, 'unique_association_fields': {'phenotype': 'http://www.ebi.ac.uk/efo/EFO_0003840', 'alleleOrigin': 'somatic', 'gene': 'ENSG00000135486', 'clinvarAccession': 'RCV000122455'}}

        test_ev_string = evidence_strings.CTTVEvidenceString(test_dict)

        self.assertEqual(self.evidence_string, test_ev_string)


class GetCTTVVariantTypeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        record_single_a = \
            ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG", "alternate": "C"},
             "snp single")
        record_single_b = ({"reference": "A", "alternate": "C"}, "snp single")
        record_single_c = ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG",
                            "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACG"}, "snp single")

        cls.test_records_singles = [record_single_a, record_single_b, record_single_c]

        record_structurals_a = \
            ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
              "alternate": "C"},
             "structural variant")
        record_structurals_b = \
            ({"reference": "A",
              "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"},
             "structural variant")
        record_structurals_c = \
            ({"reference": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
              "alternate": "AGAGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"},
             "structural variant")

        cls.test_records_structurals = \
            [record_structurals_a, record_structurals_b, record_structurals_c]

    def test_get_cttv_variant_type_singles(self):
        for record in self.test_records_singles:
            self.assertEqual(evidence_strings.get_cttv_variant_type(record[0]["reference"],
                                                      record[0]["alternate"]),
                             record[1])

    def test_get_cttv_variant_type_structurals(self):
        for record in self.test_records_structurals:
            self.assertEqual(evidence_strings.get_cttv_variant_type(record[0]["reference"],
                                                      record[0]["alternate"]),
                             record[1])


class CTTVGeneticsEvidenceStringTest(unittest.TestCase):
    def setUp(self):
        self.test_args = get_args_CTTVGeneticsEvidenceString_init()
        self.test_ges = evidence_strings.CTTVGeneticsEvidenceString(*self.test_args)

    # CTTVEvidenceString tests

    def test_unique_association_field(self):
        uaf_1 = ("gene", "test_gene")
        uaf_2 = ("clinvarAccession", "test_clinvar")
        uaf_3 = ("alleleOrigin", "germline")
        uaf_4 = ("phenotype", "test_phenotype")

        self.test_ges.add_unique_association_field(*uaf_1)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_1[0]], uaf_1[1])
        self.test_ges.add_unique_association_field(*uaf_2)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_2[0]], uaf_2[1])

        self.test_ges.add_unique_association_field(*uaf_3)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_3[0]], uaf_3[1])
        self.test_ges.add_unique_association_field(*uaf_4)
        self.assertEqual(self.test_ges['unique_association_fields'][uaf_4[0]], uaf_4[1])

    def test_set_target(self):
        target = ("http://identifiers.org/ensembl/ENSG00000135486",
                  "http://identifiers.org/cttv.activity/predicted_damaging")
        self.test_ges._clear_target()
        self.test_ges.set_target(*target)
        self.assertEqual(self.test_ges['target']['id'], [target[0]])
        self.assertEqual(self.test_ges['target']['activity'], target[1])

    def test_disease(self):
        disease_id = "Ciliary dyskinesia, primary, 26"

        self.test_ges.disease = disease_id
        self.assertEqual(self.test_ges.disease, efo_term.EFOTerm(disease_id))

    def test_evidence_codes(self):
        evidence_codes = ["http://purl.obolibrary.org/obo/ECO_0000205"]
        self.test_ges.evidence_codes = evidence_codes
        self.assertEqual(self.test_ges['evidence']['evidence_codes'], evidence_codes)
        self.assertEqual(self.test_ges.evidence_codes, evidence_codes)

    def test_top_level_literature(self):
        literature = ["http://europepmc.org/abstract/MED/20301537"]
        self.test_ges.top_level_literature = literature
        self.assertEqual(self.test_ges['literature']['references'],
                         [{"lit_id": literature_id} for literature_id in literature])
        self.assertEqual(self.test_ges.top_level_literature,
                         [{"lit_id": literature_id} for literature_id in literature])

    ###

    def test_db_xref_url(self):
        url = "http://identifiers.org/clinvar.record/RCV000128628"
        self.test_ges.db_xref_url = url
        self.assertEqual(
            self.test_ges['evidence']['gene2variant']['provenance_type']['database']['dbxref']['url'],
            url)
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['database']['dbxref']['url'],
            url)
        self.assertEqual(self.test_ges.db_xref_url, url)

    def test_url(self):
        url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
        self.test_ges.url = url
        self.assertEqual(self.test_ges['evidence']['gene2variant']['urls'][0]['url'], url)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['urls'][0]['url'], url)
        self.assertEqual(self.test_ges.url, url)

    def test_gene_2_var_ev_codes(self):
        ev_codes = ['http://identifiers.org/eco/cttv_mapping_pipeline']
        self.test_ges.gene_2_var_ev_codes = ev_codes
        self.assertEqual(self.test_ges['evidence']['gene2variant']['evidence_codes'], ev_codes)
        self.assertEqual(self.test_ges.gene_2_var_ev_codes, ev_codes)

    def test_gene_2_var_func_consequence(self):
        functional_consequence = 'http://purl.obolibrary.org/obo/SO_0001583'
        self.test_ges.gene_2_var_func_consequence = functional_consequence
        self.assertEqual(self.test_ges['evidence']['gene2variant']['functional_consequence'],
                         functional_consequence)
        self.assertEqual(self.test_ges.gene_2_var_func_consequence, functional_consequence)

    def test_set_var_2_disease_literature_a(self):
        self.test_ges['evidence']['variant2disease']['provenance_type']['literature'] = {}

        literature_1 = "PMCID12345"
        self.test_ges.set_var_2_disease_literature([literature_1])
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'],
            [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ges.set_var_2_disease_literature(literature_list)
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'],
            [{"lit_id": literature_id} for literature_id in literature_list])

    def test_set_var_2_disease_literature_b(self):
        literature_1 = "PMCID12345"
        self.test_ges.set_var_2_disease_literature([literature_1])
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'],
            [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ges.set_var_2_disease_literature(literature_list)
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['provenance_type']['literature']['references'],
            [{"lit_id": literature_id} for literature_id in literature_list])

    def test_association(self):
        self.test_ges.association = True
        self.assertTrue(self.test_ges['evidence']['gene2variant']['is_associated'])
        self.assertTrue(self.test_ges['evidence']['variant2disease']['is_associated'])
        self.assertTrue(self.test_ges.association)

        self.test_ges.association = False
        self.assertFalse(self.test_ges['evidence']['gene2variant']['is_associated'])
        self.assertFalse(self.test_ges['evidence']['variant2disease']['is_associated'])
        self.assertFalse(self.test_ges.association)

    def test_set_variant(self):
        test_id = "http://identifiers.org/dbsnp/rs193922494"
        test_type = "snp single"
        self.test_ges._clear_variant()
        self.test_ges.set_variant(test_id, test_type)
        self.assertEqual(self.test_ges['variant']['id'], [test_id])
        self.assertEqual(self.test_ges['variant']['type'], test_type)

    def test_unique_reference(self):
        unique_reference = "http://europepmc.org/abstract/MED/0"
        self.test_ges.unique_reference = unique_reference
        self.assertEqual(
            self.test_ges['evidence']['variant2disease']['unique_experiment_reference'],
            unique_reference)
        self.assertEqual(self.test_ges.unique_reference, unique_reference)

    def test_date(self):
        date_string = datetime.fromtimestamp(1412982000000 / 1000).isoformat()
        self.test_ges.date = date_string
        self.assertEqual(self.test_ges['evidence']['gene2variant']['date_asserted'], date_string)
        self.assertEqual(self.test_ges['evidence']['variant2disease']['date_asserted'],
                         date_string)
        self.assertEqual(self.test_ges.date, date_string)

    def test_validate(self):
        test_args = get_args_CTTVGeneticsEvidenceString_init()
        test_evidence_string = evidence_strings.CTTVGeneticsEvidenceString(*test_args)
        self.assertTrue(test_evidence_string.validate())


class CTTVSomaticEvidenceStringTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.consequence_type_dict = consequence_type.process_consequence_type_file(test_config.snp_2_gene_file)

    def setUp(self):
        test_args = get_args_CTTVSomaticEvidenceString_init()
        self.test_ses = evidence_strings.CTTVSomaticEvidenceString(*test_args)

    def test_db_xref_url(self):
        url = "http://identifiers.org/clinvar.record/RCV000128628"
        self.test_ses.db_xref_url = url
        self.assertEqual(self.test_ses['evidence']['provenance_type']['database']['dbxref']['url'],
                         url)
        self.assertEqual(self.test_ses.db_xref_url, url)

    def test_url(self):
        url = "http://www.ncbi.nlm.nih.gov/clinvar/RCV000128628"
        self.test_ses.url = url
        self.assertEqual(self.test_ses['evidence']['urls'][0]['url'], url)
        self.assertEqual(self.test_ses.url, url)

    def test_evidence_literature(self):
        literature_1 = "PMCID12345"
        self.test_ses.evidence_literature = [literature_1]
        self.assertEqual(self.test_ses['evidence']['provenance_type']['literature']['references'],
                         [{"lit_id": literature_1}])
        self.assertEqual(self.test_ses.evidence_literature, [{"lit_id": literature_1}])

        literature_2 = "PMCID9876"
        literature_3 = "PMCID7654"
        literature_list = [literature_2, literature_3]
        self.test_ses.evidence_literature = literature_list
        self.assertEqual(self.test_ses['evidence']['provenance_type']['literature']['references'],
                         [{"lit_id": literature_id} for literature_id in literature_list])
        self.assertEqual(self.test_ses.evidence_literature,
                         [{"lit_id": literature_id} for literature_id in literature_list])

    def test_association(self):
        self.test_ses.association = True
        self.assertTrue(self.test_ses['evidence']['is_associated'])
        self.assertTrue(self.test_ses.association)

        self.test_ses.association = False
        self.assertFalse(self.test_ses['evidence']['is_associated'])
        self.assertFalse(self.test_ses.association)

    def test_date(self):
        date_string = datetime.fromtimestamp(1412982000000 / 1000).isoformat()
        self.test_ses.date = date_string
        self.assertEqual(self.test_ses['evidence']['date_asserted'], date_string)
        self.assertEqual(self.test_ses.date, date_string)

    def test_add_known_mutations(self):
        functional_consequence = "http://purl.obolibrary.org/obo/SO_0001791"
        preferred_name = "exon_variant"
        self.test_ses._clear_known_mutations()
        self.test_ses.add_known_mutation(functional_consequence, preferred_name)
        self.assertEqual(
            self.test_ses['evidence']['known_mutations'],
            [{'functional_consequence': functional_consequence, 'preferred_name': preferred_name}])

    def test_set_known_mutations(self):
        test_consequence_type = consequence_type.ConsequenceType(ensembl_gene_ids=["ENSG00000008710"],
                                                   so_names=["3_prime_UTR_variant", "not_in_dict"])
        self.test_ses._clear_known_mutations()
        self.test_ses.set_known_mutations(test_consequence_type)
        self.assertEqual(
            self.test_ses['evidence']['known_mutations'],
            [{'functional_consequence': 'http://purl.obolibrary.org/obo/SO_0001624',
              'preferred_name': '3_prime_UTR_variant'},
             {'functional_consequence': 'http://targetvalidation.org/sequence/not_in_dict',
              'preferred_name': 'not_in_dict'}])

    def test_validate(self):
        test_args = get_args_CTTVSomaticEvidenceString_init()
        test_evidence_string = evidence_strings.CTTVSomaticEvidenceString(*test_args)
        self.assertTrue(test_evidence_string.validate())
