from datetime import datetime
import unittest
from types import SimpleNamespace

from datadiff import diff

from eva_cttv_pipeline import clinvar, efo_term, clinvar_to_evidence_strings, \
    evidence_strings, consequence_type

import tests.test_config as test_config
from tests import test_clinvar_to_evidence_strings

DATE_LAST_UPDATED = 1412982000000
DATE_ASSERTED = datetime.fromtimestamp((DATE_LAST_UPDATED / 1000)).isoformat()


MAPPINGS = SimpleNamespace()


def _get_test_cellbase_record_gene():
    return {"chromosome":"3","start":150645894,"end":150645894,"reference":"A","alternate":"C","source":"clinvar","clinvarSet":{"recordStatus":"current","title":"NM_001195794.1(CLRN1):c.567T>G (p.Tyr189Ter) AND Usher syndrome, type 3","referenceClinVarAssertion":{"clinVarAccession":{"acc":"RCV000004642","version":3,"type":"RCV","dateUpdated":1435446000000},"recordStatus":"current","clinicalSignificance":{"reviewStatus":"CRITERIA_PROVIDED_SINGLE_SUBMITTER","description":"Pathogenic","dateLastEvaluated":1435100400000},"assertion":{"type":"VARIATION_TO_DISEASE"},"attributeSet":[{"attribute":{"value":"Autosomal recessive inheritance","integerValue":263,"type":"ModeOfInheritance"},"xref":[{"db":"Laboratory for Molecular Medicine, Partners HealthCare Personalized Medicine","id":"11483565","status":"CURRENT"}]}],"observedIn":[{"sample":{"origin":"germline","species":{"value":"human","taxonomyId":9606},"affectedStatus":"not provided"},"method":[{"methodType":"CLINICAL_TESTING"},{"methodType":"LITERATURE_ONLY"}],"observedData":[{"attribute":{"integerValue":2,"type":"NumFamiliesWithVariant"},"id":6557764},{"attribute":{"value":"not provided","type":"Description"},"id":6557764},{"attribute":{"value":"Fields et al. (2002) demonstrated that the Fin(major) USH3A mutation in exon 3 of the USH3A gene, which had been identified by Joensuu et al. (2001) as 300C-T (TYR100TER), should be referred to as 528T-G, resulting in a tyr176-to-ter substitution. Joensuu et al. (2001) had identified homozygosity for this mutation in a Finnish family segregating Usher syndrome type IIIA (USH3A; 276902) and found it in a further 52 Finnish patients. Fields et al. (2002) found this mutation in 11 of 28 mutated alleles from affected individuals of Finnish and other northern European ancestry.","type":"Description"},"citation":[{"id":[{"value":"11524702","source":"PubMed"}],"type":"general"},{"id":[{"value":"12145752","source":"PubMed"}],"type":"general"}],"id":6557764},{"attribute":{"integerValue":3,"type":"VariantAlleles"},"id":6557764}]}],"measureSet":{"measure":[{"name":[{"elementValue":{"value":"NM_001195794.1(CLRN1):c.567T>G (p.Tyr189Ter)","type":"Preferred"}}],"attributeSet":[{"attribute":{"value":"0.000076887590","type":"AlleleFrequency"},"xref":[{"db":"dbSNP","id":"121908140","status":"CURRENT"},{"db":"NHLBI GO Exome Sequencing Project (ESP)","id":"ESP6500SI-V2","url":"http://evs.gs.washington.edu/EVS/","status":"CURRENT"}]},{"attribute":{"value":"NM_174878.2:c.528T>G","type":"HGVS, coding","change":"c.528T>G"}},{"attribute":{"value":"NM_001256819.1:c.*142T>G","type":"HGVS, coding, RefSeq","change":"c.*142T>G"}},{"attribute":{"value":"NM_052995.2:c.300T>G","type":"HGVS, coding, RefSeq","change":"c.300T>G"}},{"attribute":{"value":"NM_001195794.1:c.567T>G","type":"HGVS, coding, RefSeq","change":"c.567T>G"}},{"attribute":{"value":"NG_009168.1:g.49893T>G","type":"HGVS, genomic, RefSeqGene","change":"g.49893T>G"}},{"attribute":{"value":"NC_000003.12:g.150928107A>C","integerValue":38,"type":"HGVS, genomic, top level","change":"g.150928107A>C"}},{"attribute":{"value":"NC_000003.11:g.150645894A>C","integerValue":37,"type":"HGVS, genomic, top level, previous","change":"g.150645894A>C"}},{"attribute":{"value":"NR_046380.2:n.1009T>G","type":"HGVS, non-coding","change":"n.1009T>G"}},{"attribute":{"value":"NR_046380.1:n.1010T>G","type":"HGVS, previous","change":"n.1010T>G"}},{"attribute":{"value":"p.Tyr176X","type":"HGVS, protein"}},{"attribute":{"value":"NP_443721.1:p.Tyr100Ter","type":"HGVS, protein, RefSeq","change":"p.Tyr100Ter"},"xref":[{"db":"dbSNP","id":"121908140","type":"rs","status":"CURRENT"}]},{"attribute":{"value":"NP_777367.1:p.Tyr176Ter","type":"HGVS, protein, RefSeq","change":"p.Tyr176Ter"},"xref":[{"db":"dbSNP","id":"121908140","type":"rs","status":"CURRENT"}]},{"attribute":{"value":"NP_001182723.1:p.Tyr189Ter","type":"HGVS, protein, RefSeq","change":"p.Tyr189Ter"},"xref":[{"db":"dbSNP","id":"121908140","type":"rs","status":"CURRENT"}]},{"attribute":{"value":"NM_174878.2:EXON 3","type":"Location"}},{"attribute":{"value":"3 prime UTR variant","type":"MolecularConsequence"},"xref":[{"db":"Sequence Ontology","id":"SO:0001624","status":"CURRENT"},{"db":"RefSeq","id":"NM_001256819.1:c.*142T>G","status":"CURRENT"}]},{"attribute":{"value":"nonsense","type":"MolecularConsequence"},"xref":[{"db":"Sequence Ontology","id":"SO:0001587","status":"CURRENT"},{"db":"RefSeq","id":"NM_001195794.1:c.567T>G","status":"CURRENT"}]},{"attribute":{"value":"non-coding transcript variant","type":"MolecularConsequence"},"xref":[{"db":"Sequence Ontology","id":"SO:0001619","status":"CURRENT"},{"db":"RefSeq","id":"NR_046380.2:n.1009T>G","status":"CURRENT"}]},{"attribute":{"value":"Y176*","type":"ProteinChange1LetterCode"},"xref":[{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"}]},{"attribute":{"value":"Y100*","type":"ProteinChange1LetterCode"}},{"attribute":{"value":"Y189*","type":"ProteinChange1LetterCode"}},{"attribute":{"value":"TYR176TER","type":"ProteinChange3LetterCode"},"xref":[{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"}]}],"cytogeneticLocation":["3q25.1"],"sequenceLocation":[{"assembly":"GRCh38","chr":"3","accession":"NC_000003.12","start":150928107,"stop":150928107,"displayStart":150928107,"displayStop":150928107,"variantLength":1,"referenceAllele":"A","alternateAllele":"C","assemblyAccessionVersion":"GCF_000001405.26","assemblyStatus":"current"},{"assembly":"GRCh37","chr":"3","accession":"NC_000003.11","start":150645894,"stop":150645894,"displayStart":150645894,"displayStop":150645894,"variantLength":1,"referenceAllele":"A","alternateAllele":"C","assemblyAccessionVersion":"GCF_000001405.25","assemblyStatus":"previous"}],"measureRelationship":[{"name":[{"elementValue":{"value":"clarin 1","type":"Preferred"}}],"symbol":[{"elementValue":{"value":"CLRN1","type":"Preferred"}}],"sequenceLocation":[{"assembly":"GRCh38","chr":"3","accession":"NC_000003.12","start":150918910,"stop":150973019,"displayStart":150918910,"displayStop":150973019,"strand":"-","variantLength":46837,"assemblyAccessionVersion":"GCF_000001405.26","assemblyStatus":"current"},{"assembly":"GRCh37","chr":"3","accession":"NC_000003.11","start":150643949,"stop":150690785,"displayStart":150643949,"displayStop":150690785,"strand":"-","variantLength":46837,"assemblyAccessionVersion":"GCF_000001405.25","assemblyStatus":"previous"}],"type":"variant in gene","xref":[{"db":"Gene","id":"7401","status":"CURRENT"},{"db":"OMIM","id":"606397","type":"MIM","status":"CURRENT"}]}],"type":"single nucleotide variant","id":19431,"xref":[{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"},{"db":"dbSNP","id":"121908140","type":"rs","status":"CURRENT"}]}],"name":[{"elementValue":{"value":"NM_001195794.1(CLRN1):c.567T>G (p.Tyr189Ter)","type":"Preferred"}}],"type":"Variant","id":4392},"traitSet":{"trait":[{"name":[{"elementValue":{"value":"Usher syndrome, type 3","type":"Preferred"},"xref":[{"db":"Genetic Alliance","id":"Usher+syndrome%2C+type+3/7326","status":"CURRENT"},{"db":"Office of Rare Diseases","id":"5442","status":"CURRENT"}]},{"elementValue":{"value":"Usher Syndrome, Type III","type":"Alternate"}},{"elementValue":{"value":"USHER SYNDROME, TYPE IIIA","type":"Alternate"},"xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"},{"db":"OMIM","id":"606397.0002","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0007","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0004","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0005","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0003","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0008","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0006","type":"Allelic variant","status":"CURRENT"}]},{"elementValue":{"value":"Usher syndrome, type 3A","type":"Alternate"}},{"elementValue":{"value":"Orphanet:886","type":"EFO id"}},{"elementValue":{"value":"Usher syndrome","type":"EFO name"}},{"elementValue":{"value":"http://www.orpha.net/ORDO/Orphanet_886","type":"EFO URL"}}],"symbol":[{"elementValue":{"value":"USH3","type":"Preferred"},"xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"},{"db":"Office of Rare Diseases","id":"5442","status":"CURRENT"}]},{"elementValue":{"value":"USH3A","type":"Alternate"},"xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"},{"db":"Office of Rare Diseases","id":"5442","status":"CURRENT"}]}],"attributeSet":[{"attribute":{"value":"Neonatal/infancy","type":"age of onset"},"xref":[{"db":"Orphanet","id":"886","status":"CURRENT"},{"db":"Orphanet","id":"231183","status":"CURRENT"}]}],"citation":[{"id":[{"value":"21697857","source":"PubMed"}],"type":"Translational/Evidence-based","abbrev":"EuroGenetest, 2011"}],"type":"Disease","id":5092,"xref":[{"db":"MedGen","id":"C1568248","status":"CURRENT"},{"db":"Orphanet","id":"231183","status":"CURRENT"},{"db":"Orphanet","id":"886","status":"CURRENT"},{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"}]}],"type":"Disease","id":1209},"dateCreated":1344812400000,"dateLastUpdated":1435359600000,"id":62145},"clinVarAssertion":[{"clinVarSubmissionID":{"submitter":"OMIM","title":"CLRN1, TYR176TER_USHER SYNDROME, TYPE IIIA","localKey":"606397.0001_USHER SYNDROME, TYPE IIIA","submitterDate":1435100400000},"clinVarAccession":{"acc":"SCV000024816","version":2,"type":"SCV","orgID":3,"dateUpdated":1435359600000},"recordStatus":"current","clinicalSignificance":{"reviewStatus":"NO_ASSERTION_CRITERIA_PROVIDED","description":["Pathogenic"],"dateLastEvaluated":1435100400000},"assertion":{"type":"variation to disease"},"externalID":{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"},"observedIn":[{"sample":{"origin":"germline","species":{"value":"human"},"affectedStatus":"not provided"},"method":[{"methodType":"LITERATURE_ONLY"}],"observedData":[{"attribute":{"value":"Fields et al. (2002) demonstrated that the Fin(major) USH3A mutation in exon 3 of the USH3A gene, which had been identified by Joensuu et al. (2001) as 300C-T (TYR100TER), should be referred to as 528T-G, resulting in a tyr176-to-ter substitution. Joensuu et al. (2001) had identified homozygosity for this mutation in a Finnish family segregating Usher syndrome type IIIA (USH3A; 276902) and found it in a further 52 Finnish patients. Fields et al. (2002) found this mutation in 11 of 28 mutated alleles from affected individuals of Finnish and other northern European ancestry.","type":"Description"},"citation":[{"id":[{"value":"12145752","source":"PubMed"}]},{"id":[{"value":"11524702","source":"PubMed"}]}],"xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"}]}]}],"measureSet":{"measure":[{"name":[{"elementValue":{"value":"CLRN1, TYR176TER","type":"Preferred"}}],"attributeSet":[{"attribute":{"value":"TYR176TER","type":"NonHGVS"}}],"measureRelationship":[{"symbol":[{"elementValue":{"value":"CLRN1","type":"Preferred"}}],"type":"variant in gene"}],"type":"Variation","xref":[{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"}]}],"type":"Variant"},"traitSet":{"trait":[{"name":[{"elementValue":{"value":"USHER SYNDROME, TYPE IIIA","type":"Preferred"}}],"type":"Disease"}],"type":"Disease"},"id":24816},{"clinVarSubmissionID":{"submitter":"Laboratory for Molecular Medicine,Partners HealthCare Personalized Medicine","localKey":"11483565|OMIM:276902","submitterDate":1422489600000},"clinVarAccession":{"acc":"SCV000203992","version":1,"type":"SCV","orgID":21766,"dateUpdated":1422576000000},"recordStatus":"current","clinicalSignificance":{"reviewStatus":"CLASSIFIED_BY_SINGLE_SUBMITTER","description":["Pathogenic"],"citation":[{"id":[{"value":"11524702","source":"PubMed"}]}],"comment":[{"value":"The Tyr176X variant in CLRN1 has been previously identified in 52 homozygous and 2 compound heterozygous individuals with Usher syndrome type III (Joensuu 2001). This variant has been identified in 1/8,600 European American chromosomes by the NHLBI Exome Sequencing Project (http://evs.gs.washington.edu/EVS/; dbSNP rs121908140). Although this variant has been seen in the general population, its frequency is low enough to be consistent with a recessive carrier frequency. This nonsense variant leads to a premature termination codon at position 176, which is predicted to lead to a truncated or absent protein. In summary, this variant meets our criteria to be classified as pathogenic in a recessive manner for Usher syndrome (http://pcpgm.partners.org/LMM)."}],"dateLastEvaluated":1388620800000},"assertion":{"type":"variation to disease"},"externalID":{"db":"Laboratory for Molecular Medicine (Partners HealthCare Personalized Medicine)","id":"11483565","status":"CURRENT"},"attributeSet":[{"attribute":{"value":"Autosomal recessive inheritance","type":"ModeOfInheritance"}}],"observedIn":[{"sample":{"origin":"germline","species":{"value":"human","taxonomyId":9606},"affectedStatus":"not provided","familyData":{"numFamiliesWithVariant":2}},"method":[{"methodType":"CLINICAL_TESTING"}],"observedData":[{"attribute":{"integerValue":3,"type":"VariantAlleles"}}]}],"measureSet":{"measure":[{"name":[{"elementValue":{"value":"NM_174878.2:c.528T>G","type":"Alternate"}},{"elementValue":{"value":"p.Tyr176X","type":"Alternate"}}],"attributeSet":[{"attribute":{"value":"NM_174878.2:EXON 3","type":"Location"}},{"attribute":{"value":"NC_000003.11:g.150645894A>C","type":"HGVS"}}],"sequenceLocation":[{"assembly":"GRCh37","chr":"3","start":150645894,"stop":150645894,"variantLength":1,"referenceAllele":"A","alternateAllele":"C"}],"measureRelationship":[{"symbol":[{"elementValue":{"value":"CLRN1","type":"Preferred"}}],"type":"variant in gene"}],"type":"Variation","xref":[{"db":"dbSNP","id":"121908140","type":"rsNumber","status":"CURRENT"}]}],"type":"Variant"},"traitSet":{"trait":[{"name":[{"elementValue":{"value":"Usher syndrome, type 3A","type":"Preferred"}}],"type":"Disease","xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"}]}],"type":"Disease"},"submissionName":"LMM_all.variants_NCBI_3.16.2013","id":366075}],"id":6973966}}



def get_args_CTTVGeneticsEvidenceString_init():
    cellbase_record = _get_test_cellbase_record_gene()

    clinvarRecord = \
        clinvar.ClinvarRecord(mappings=test_clinvar_to_evidence_strings.MAPPINGS,
                                     a_dictionary={"recordStatus":"current","title":"NM_001195794.1(CLRN1):c.567T>G (p.Tyr189Ter) AND Usher syndrome, type 3","referenceClinVarAssertion":{"clinVarAccession":{"acc":"RCV000004642","version":3,"type":"RCV","dateUpdated":1435446000000},"recordStatus":"current","clinicalSignificance":{"reviewStatus":"CRITERIA_PROVIDED_SINGLE_SUBMITTER","description":"Pathogenic","dateLastEvaluated":1435100400000},"assertion":{"type":"VARIATION_TO_DISEASE"},"attributeSet":[{"attribute":{"value":"Autosomal recessive inheritance","integerValue":263,"type":"ModeOfInheritance"},"xref":[{"db":"Laboratory for Molecular Medicine, Partners HealthCare Personalized Medicine","id":"11483565","status":"CURRENT"}]}],"observedIn":[{"sample":{"origin":"germline","species":{"value":"human","taxonomyId":9606},"affectedStatus":"not provided"},"method":[{"methodType":"CLINICAL_TESTING"},{"methodType":"LITERATURE_ONLY"}],"observedData":[{"attribute":{"integerValue":2,"type":"NumFamiliesWithVariant"},"id":6557764},{"attribute":{"value":"not provided","type":"Description"},"id":6557764},{"attribute":{"value":"Fields et al. (2002) demonstrated that the Fin(major) USH3A mutation in exon 3 of the USH3A gene, which had been identified by Joensuu et al. (2001) as 300C-T (TYR100TER), should be referred to as 528T-G, resulting in a tyr176-to-ter substitution. Joensuu et al. (2001) had identified homozygosity for this mutation in a Finnish family segregating Usher syndrome type IIIA (USH3A; 276902) and found it in a further 52 Finnish patients. Fields et al. (2002) found this mutation in 11 of 28 mutated alleles from affected individuals of Finnish and other northern European ancestry.","type":"Description"},"citation":[{"id":[{"value":"11524702","source":"PubMed"}],"type":"general"},{"id":[{"value":"12145752","source":"PubMed"}],"type":"general"}],"id":6557764},{"attribute":{"integerValue":3,"type":"VariantAlleles"},"id":6557764}]}],"measureSet":{"measure":[{"name":[{"elementValue":{"value":"NM_001195794.1(CLRN1):c.567T>G (p.Tyr189Ter)","type":"Preferred"}}],"attributeSet":[{"attribute":{"value":"0.000076887590","type":"AlleleFrequency"},"xref":[{"db":"dbSNP","id":"121908140","status":"CURRENT"},{"db":"NHLBI GO Exome Sequencing Project (ESP)","id":"ESP6500SI-V2","url":"http://evs.gs.washington.edu/EVS/","status":"CURRENT"}]},{"attribute":{"value":"NM_174878.2:c.528T>G","type":"HGVS, coding","change":"c.528T>G"}},{"attribute":{"value":"NM_001256819.1:c.*142T>G","type":"HGVS, coding, RefSeq","change":"c.*142T>G"}},{"attribute":{"value":"NM_052995.2:c.300T>G","type":"HGVS, coding, RefSeq","change":"c.300T>G"}},{"attribute":{"value":"NM_001195794.1:c.567T>G","type":"HGVS, coding, RefSeq","change":"c.567T>G"}},{"attribute":{"value":"NG_009168.1:g.49893T>G","type":"HGVS, genomic, RefSeqGene","change":"g.49893T>G"}},{"attribute":{"value":"NC_000003.12:g.150928107A>C","integerValue":38,"type":"HGVS, genomic, top level","change":"g.150928107A>C"}},{"attribute":{"value":"NC_000003.11:g.150645894A>C","integerValue":37,"type":"HGVS, genomic, top level, previous","change":"g.150645894A>C"}},{"attribute":{"value":"NR_046380.2:n.1009T>G","type":"HGVS, non-coding","change":"n.1009T>G"}},{"attribute":{"value":"NR_046380.1:n.1010T>G","type":"HGVS, previous","change":"n.1010T>G"}},{"attribute":{"value":"p.Tyr176X","type":"HGVS, protein"}},{"attribute":{"value":"NP_443721.1:p.Tyr100Ter","type":"HGVS, protein, RefSeq","change":"p.Tyr100Ter"},"xref":[{"db":"dbSNP","id":"121908140","type":"rs","status":"CURRENT"}]},{"attribute":{"value":"NP_777367.1:p.Tyr176Ter","type":"HGVS, protein, RefSeq","change":"p.Tyr176Ter"},"xref":[{"db":"dbSNP","id":"121908140","type":"rs","status":"CURRENT"}]},{"attribute":{"value":"NP_001182723.1:p.Tyr189Ter","type":"HGVS, protein, RefSeq","change":"p.Tyr189Ter"},"xref":[{"db":"dbSNP","id":"121908140","type":"rs","status":"CURRENT"}]},{"attribute":{"value":"NM_174878.2:EXON 3","type":"Location"}},{"attribute":{"value":"3 prime UTR variant","type":"MolecularConsequence"},"xref":[{"db":"Sequence Ontology","id":"SO:0001624","status":"CURRENT"},{"db":"RefSeq","id":"NM_001256819.1:c.*142T>G","status":"CURRENT"}]},{"attribute":{"value":"nonsense","type":"MolecularConsequence"},"xref":[{"db":"Sequence Ontology","id":"SO:0001587","status":"CURRENT"},{"db":"RefSeq","id":"NM_001195794.1:c.567T>G","status":"CURRENT"}]},{"attribute":{"value":"non-coding transcript variant","type":"MolecularConsequence"},"xref":[{"db":"Sequence Ontology","id":"SO:0001619","status":"CURRENT"},{"db":"RefSeq","id":"NR_046380.2:n.1009T>G","status":"CURRENT"}]},{"attribute":{"value":"Y176*","type":"ProteinChange1LetterCode"},"xref":[{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"}]},{"attribute":{"value":"Y100*","type":"ProteinChange1LetterCode"}},{"attribute":{"value":"Y189*","type":"ProteinChange1LetterCode"}},{"attribute":{"value":"TYR176TER","type":"ProteinChange3LetterCode"},"xref":[{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"}]}],"cytogeneticLocation":["3q25.1"],"sequenceLocation":[{"assembly":"GRCh38","chr":"3","accession":"NC_000003.12","start":150928107,"stop":150928107,"displayStart":150928107,"displayStop":150928107,"variantLength":1,"referenceAllele":"A","alternateAllele":"C","assemblyAccessionVersion":"GCF_000001405.26","assemblyStatus":"current"},{"assembly":"GRCh37","chr":"3","accession":"NC_000003.11","start":150645894,"stop":150645894,"displayStart":150645894,"displayStop":150645894,"variantLength":1,"referenceAllele":"A","alternateAllele":"C","assemblyAccessionVersion":"GCF_000001405.25","assemblyStatus":"previous"}],"measureRelationship":[{"name":[{"elementValue":{"value":"clarin 1","type":"Preferred"}}],"symbol":[{"elementValue":{"value":"CLRN1","type":"Preferred"}}],"sequenceLocation":[{"assembly":"GRCh38","chr":"3","accession":"NC_000003.12","start":150918910,"stop":150973019,"displayStart":150918910,"displayStop":150973019,"strand":"-","variantLength":46837,"assemblyAccessionVersion":"GCF_000001405.26","assemblyStatus":"current"},{"assembly":"GRCh37","chr":"3","accession":"NC_000003.11","start":150643949,"stop":150690785,"displayStart":150643949,"displayStop":150690785,"strand":"-","variantLength":46837,"assemblyAccessionVersion":"GCF_000001405.25","assemblyStatus":"previous"}],"type":"variant in gene","xref":[{"db":"Gene","id":"7401","status":"CURRENT"},{"db":"OMIM","id":"606397","type":"MIM","status":"CURRENT"}]}],"type":"single nucleotide variant","id":19431,"xref":[{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"},{"db":"dbSNP","id":"121908140","type":"rs","status":"CURRENT"}]}],"name":[{"elementValue":{"value":"NM_001195794.1(CLRN1):c.567T>G (p.Tyr189Ter)","type":"Preferred"}}],"type":"Variant","id":4392},"traitSet":{"trait":[{"name":[{"elementValue":{"value":"Usher syndrome, type 3","type":"Preferred"},"xref":[{"db":"Genetic Alliance","id":"Usher+syndrome%2C+type+3/7326","status":"CURRENT"},{"db":"Office of Rare Diseases","id":"5442","status":"CURRENT"}]},{"elementValue":{"value":"Usher Syndrome, Type III","type":"Alternate"}},{"elementValue":{"value":"USHER SYNDROME, TYPE IIIA","type":"Alternate"},"xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"},{"db":"OMIM","id":"606397.0002","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0007","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0004","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0005","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0003","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0008","type":"Allelic variant","status":"CURRENT"},{"db":"OMIM","id":"606397.0006","type":"Allelic variant","status":"CURRENT"}]},{"elementValue":{"value":"Usher syndrome, type 3A","type":"Alternate"}},{"elementValue":{"value":"Orphanet:886","type":"EFO id"}},{"elementValue":{"value":"Usher syndrome","type":"EFO name"}},{"elementValue":{"value":"http://www.orpha.net/ORDO/Orphanet_886","type":"EFO URL"}}],"symbol":[{"elementValue":{"value":"USH3","type":"Preferred"},"xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"},{"db":"Office of Rare Diseases","id":"5442","status":"CURRENT"}]},{"elementValue":{"value":"USH3A","type":"Alternate"},"xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"},{"db":"Office of Rare Diseases","id":"5442","status":"CURRENT"}]}],"attributeSet":[{"attribute":{"value":"Neonatal/infancy","type":"age of onset"},"xref":[{"db":"Orphanet","id":"886","status":"CURRENT"},{"db":"Orphanet","id":"231183","status":"CURRENT"}]}],"citation":[{"id":[{"value":"21697857","source":"PubMed"}],"type":"Translational/Evidence-based","abbrev":"EuroGenetest, 2011"}],"type":"Disease","id":5092,"xref":[{"db":"MedGen","id":"C1568248","status":"CURRENT"},{"db":"Orphanet","id":"231183","status":"CURRENT"},{"db":"Orphanet","id":"886","status":"CURRENT"},{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"}]}],"type":"Disease","id":1209},"dateCreated":1344812400000,"dateLastUpdated":1435359600000,"id":62145},"clinVarAssertion":[{"clinVarSubmissionID":{"submitter":"OMIM","title":"CLRN1, TYR176TER_USHER SYNDROME, TYPE IIIA","localKey":"606397.0001_USHER SYNDROME, TYPE IIIA","submitterDate":1435100400000},"clinVarAccession":{"acc":"SCV000024816","version":2,"type":"SCV","orgID":3,"dateUpdated":1435359600000},"recordStatus":"current","clinicalSignificance":{"reviewStatus":"NO_ASSERTION_CRITERIA_PROVIDED","description":["Pathogenic"],"dateLastEvaluated":1435100400000},"assertion":{"type":"variation to disease"},"externalID":{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"},"observedIn":[{"sample":{"origin":"germline","species":{"value":"human"},"affectedStatus":"not provided"},"method":[{"methodType":"LITERATURE_ONLY"}],"observedData":[{"attribute":{"value":"Fields et al. (2002) demonstrated that the Fin(major) USH3A mutation in exon 3 of the USH3A gene, which had been identified by Joensuu et al. (2001) as 300C-T (TYR100TER), should be referred to as 528T-G, resulting in a tyr176-to-ter substitution. Joensuu et al. (2001) had identified homozygosity for this mutation in a Finnish family segregating Usher syndrome type IIIA (USH3A; 276902) and found it in a further 52 Finnish patients. Fields et al. (2002) found this mutation in 11 of 28 mutated alleles from affected individuals of Finnish and other northern European ancestry.","type":"Description"},"citation":[{"id":[{"value":"12145752","source":"PubMed"}]},{"id":[{"value":"11524702","source":"PubMed"}]}],"xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"}]}]}],"measureSet":{"measure":[{"name":[{"elementValue":{"value":"CLRN1, TYR176TER","type":"Preferred"}}],"attributeSet":[{"attribute":{"value":"TYR176TER","type":"NonHGVS"}}],"measureRelationship":[{"symbol":[{"elementValue":{"value":"CLRN1","type":"Preferred"}}],"type":"variant in gene"}],"type":"Variation","xref":[{"db":"OMIM","id":"606397.0001","type":"Allelic variant","status":"CURRENT"}]}],"type":"Variant"},"traitSet":{"trait":[{"name":[{"elementValue":{"value":"USHER SYNDROME, TYPE IIIA","type":"Preferred"}}],"type":"Disease"}],"type":"Disease"},"id":24816},{"clinVarSubmissionID":{"submitter":"Laboratory for Molecular Medicine,Partners HealthCare Personalized Medicine","localKey":"11483565|OMIM:276902","submitterDate":1422489600000},"clinVarAccession":{"acc":"SCV000203992","version":1,"type":"SCV","orgID":21766,"dateUpdated":1422576000000},"recordStatus":"current","clinicalSignificance":{"reviewStatus":"CLASSIFIED_BY_SINGLE_SUBMITTER","description":["Pathogenic"],"citation":[{"id":[{"value":"11524702","source":"PubMed"}]}],"comment":[{"value":"The Tyr176X variant in CLRN1 has been previously identified in 52 homozygous and 2 compound heterozygous individuals with Usher syndrome type III (Joensuu 2001). This variant has been identified in 1/8,600 European American chromosomes by the NHLBI Exome Sequencing Project (http://evs.gs.washington.edu/EVS/; dbSNP rs121908140). Although this variant has been seen in the general population, its frequency is low enough to be consistent with a recessive carrier frequency. This nonsense variant leads to a premature termination codon at position 176, which is predicted to lead to a truncated or absent protein. In summary, this variant meets our criteria to be classified as pathogenic in a recessive manner for Usher syndrome (http://pcpgm.partners.org/LMM)."}],"dateLastEvaluated":1388620800000},"assertion":{"type":"variation to disease"},"externalID":{"db":"Laboratory for Molecular Medicine (Partners HealthCare Personalized Medicine)","id":"11483565","status":"CURRENT"},"attributeSet":[{"attribute":{"value":"Autosomal recessive inheritance","type":"ModeOfInheritance"}}],"observedIn":[{"sample":{"origin":"germline","species":{"value":"human","taxonomyId":9606},"affectedStatus":"not provided","familyData":{"numFamiliesWithVariant":2}},"method":[{"methodType":"CLINICAL_TESTING"}],"observedData":[{"attribute":{"integerValue":3,"type":"VariantAlleles"}}]}],"measureSet":{"measure":[{"name":[{"elementValue":{"value":"NM_174878.2:c.528T>G","type":"Alternate"}},{"elementValue":{"value":"p.Tyr176X","type":"Alternate"}}],"attributeSet":[{"attribute":{"value":"NM_174878.2:EXON 3","type":"Location"}},{"attribute":{"value":"NC_000003.11:g.150645894A>C","type":"HGVS"}}],"sequenceLocation":[{"assembly":"GRCh37","chr":"3","start":150645894,"stop":150645894,"variantLength":1,"referenceAllele":"A","alternateAllele":"C"}],"measureRelationship":[{"symbol":[{"elementValue":{"value":"CLRN1","type":"Preferred"}}],"type":"variant in gene"}],"type":"Variation","xref":[{"db":"dbSNP","id":"121908140","type":"rsNumber","status":"CURRENT"}]}],"type":"Variant"},"traitSet":{"trait":[{"name":[{"elementValue":{"value":"Usher syndrome, type 3A","type":"Preferred"}}],"type":"Disease","xref":[{"db":"OMIM","id":"276902","type":"MIM","status":"CURRENT"}]}],"type":"Disease"},"submissionName":"LMM_all.variants_NCBI_3.16.2013","id":366075}],"id":6973966})

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
        test_dict = {"literature": {"references": [{"lit_id": "http://europepmc.org/abstract/MED/12145752"}, {"lit_id": "http://europepmc.org/abstract/MED/21697857"}, {"lit_id": "http://europepmc.org/abstract/MED/11524702"}]}, "disease": {"id": ["http://www.orpha.net/ORDO/Orphanet_886"]}, "validated_against_schema_version": "1.2.3", "target": {"target_type": "http://identifiers.org/cttv.target/gene_variant", "id": ["http://identifiers.org/ensembl/ENSG00000163646"], "activity": "http://identifiers.org/cttv.activity/damaging_to_target"}, "sourceID": "eva", "evidence": {"gene2variant": {"is_associated": True, "provenance_type": {"expert": {"status": True, "statement": "Primary submitter of data"}, "database": {"id": "EVA", "dbxref": {"url": "http://identifiers.org/clinvar.record/RCV000004642", "id": "http://identifiers.org/clinvar", "version": "2015-04"}, "version": "1.0"}}, "evidence_codes": ["http://identifiers.org/eco/cttv_mapping_pipeline"], "date_asserted": "2015-06-27T00:00:00", "functional_consequence": "http://purl.obolibrary.org/obo/SO_0001587", "urls": [{"url": "http://www.ncbi.nlm.nih.gov/clinvar/RCV000004642", "nice_name": "Further details in ClinVar database"}]}, "variant2disease": {"is_associated": True, "provenance_type": {"literature": {"references": [{"lit_id": "http://europepmc.org/abstract/MED/12145752"}, {"lit_id": "http://europepmc.org/abstract/MED/21697857"}, {"lit_id": "http://europepmc.org/abstract/MED/11524702"}]}, "expert": {"status": True, "statement": "Primary submitter of data"}, "database": {"id": "EVA", "dbxref": {"url": "http://identifiers.org/clinvar.record/RCV000004642", "id": "http://identifiers.org/clinvar", "version": "2015-04"}, "version": "1.0"}}, "evidence_codes": ["http://purl.obolibrary.org/obo/ECO_0000205"], "date_asserted": "2015-06-27T00:00:00", "unique_experiment_reference": "http://europepmc.org/abstract/MED/12145752", "urls": [{"url": "http://www.ncbi.nlm.nih.gov/clinvar/RCV000004642", "nice_name": "Further details in ClinVar database"}], "resource_score": {"type": "pvalue", "method": {"url": "", "description": "Not provided by data supplier"}, "value": 1e-07}}}, "type": "genetic_association", "access_level": "public", "unique_association_fields": {"gene": "ENSG00000163646", "alleleOrigin": "germline", "phenotype": "http://www.orpha.net/ORDO/Orphanet_886", "clinvarAccession": "RCV000004642"}, "variant": {"type": "snp single", "id": ["http://identifiers.org/dbsnp/rs121908140"]}}

        test_ev_string = evidence_strings.CTTVEvidenceString(test_dict)

        self.assertEqual(set(self.evidence_string.__dict__.values()), set(test_ev_string.__dict__.values()))


def get_args_CTTVSomaticEvidenceString_init():
    clinvarRecord = clinvar.ClinvarRecord(mappings=test_clinvar_to_evidence_strings.MAPPINGS,
                                                 a_dictionary={"recordStatus":"current","title":"NM_000038.5(APC):c.4391_4394delAGAG (p.Glu1464Valfs) AND Periampullary adenoma","referenceClinVarAssertion":{"clinVarAccession":{"acc":"RCV000000851","version":4,"type":"RCV","dateUpdated":1455667200000},"recordStatus":"current","clinicalSignificance":{"reviewStatus":"NO_ASSERTION_CRITERIA_PROVIDED","description":"Pathogenic","dateLastEvaluated":752112000000},"assertion":{"type":"VARIATION_TO_DISEASE"},"observedIn":[{"sample":{"origin":"somatic","species":{"value":"human","taxonomyId":9606},"affectedStatus":"not provided"},"method":[{"methodType":"LITERATURE_ONLY"}],"observedData":[{"attribute":{"value":"In tumor tissue of a periampullary adenoma from a patient with FAP (175100), Bapat et al. (1993) identified a somatic 4-bp deletion (AGAG) at codon 1464 of the APC gene. The patient had a germline APC mutation (611731.0023).","type":"Description"},"citation":[{"id":[{"value":"8281160","source":"PubMed"}],"type":"general"}],"id":9728450}]}],"measureSet":{"measure":[{"name":[{"elementValue":{"value":"NM_000038.5(APC):c.4391_4394delAGAG (p.Glu1464Valfs)","type":"Preferred"}}],"attributeSet":[{"attribute":{"value":"NM_001127511.2:c.4337_4340delAGAG","type":"HGVS, coding, RefSeq","change":"c.4337_4340delAGAG"}},{"attribute":{"value":"NM_000038.5:c.4391_4394delAGAG","type":"HGVS, coding, RefSeq","change":"c.4391_4394delAGAG"}},{"attribute":{"value":"NM_001127510.2:c.4391_4394delAGAG","type":"HGVS, coding, RefSeq","change":"c.4391_4394delAGAG"}},{"attribute":{"value":"LRG_130:g.152465_152468delAGAG","type":"HGVS, genomic, LRG","change":"g.152465_152468delAGAG"}},{"attribute":{"value":"NG_008481.4:g.152465_152468delAGAG","type":"HGVS, genomic, RefSeqGene","change":"g.152465_152468delAGAG"},"citation":[{"id":[{"value":"8281160","source":"PubMed"}],"type":"general"}]},{"attribute":{"value":"NC_000005.10:g.112839985_112839988delAGAG","integerValue":38,"type":"HGVS, genomic, top level","change":"g.112839985_112839988delAGAG"}},{"attribute":{"value":"NC_000005.9:g.112175682_112175685delAGAG","integerValue":37,"type":"HGVS, genomic, top level, previous","change":"g.112175682_112175685delAGAG"}},{"attribute":{"value":"NM_000038.4:c.4391_4394delAGAG","type":"HGVS, previous","change":"c.4391_4394delAGAG"}},{"attribute":{"value":"LRG_130p1:p.Glu1464Valfs","type":"HGVS, protein","change":"p.Glu1464Valfs"}},{"attribute":{"value":"LRG_130p2:p.Glu1464Valfs","type":"HGVS, protein","change":"p.Glu1464Valfs"}},{"attribute":{"value":"NP_001120983.2:p.Glu1446Valfs","type":"HGVS, protein, RefSeq","change":"p.Glu1446Valfs"}},{"attribute":{"value":"NP_000029.2:p.Glu1464Valfs","type":"HGVS, protein, RefSeq","change":"p.Glu1464Valfs"}},{"attribute":{"value":"NP_001120982.1:p.Glu1464Valfs","type":"HGVS, protein, RefSeq","change":"p.Glu1464Valfs"}},{"attribute":{"value":"frameshift variant","type":"MolecularConsequence"},"xref":[{"db":"Sequence Ontology","id":"SO:0001589","status":"CURRENT"},{"db":"RefSeq","id":"NM_000038.5:c.4391_4394delAGAG","status":"CURRENT"}]}],"cytogeneticLocation":["5q22.2"],"sequenceLocation":[{"assembly":"GRCh38","chr":"5","accession":"NC_000005.10","start":112839985,"stop":112839988,"displayStart":112839985,"displayStop":112839988,"variantLength":4,"referenceAllele":"AGAG","alternateAllele":"-","assemblyAccessionVersion":"GCF_000001405.28","assemblyStatus":"current"},{"assembly":"GRCh37","chr":"5","accession":"NC_000005.9","start":112175682,"stop":112175685,"displayStart":112175682,"displayStop":112175685,"variantLength":4,"referenceAllele":"AGAG","alternateAllele":"-","assemblyAccessionVersion":"GCF_000001405.25","assemblyStatus":"previous"}],"measureRelationship":[{"name":[{"elementValue":{"value":"adenomatous polyposis coli","type":"Preferred"}}],"symbol":[{"elementValue":{"value":"APC","type":"Preferred"}}],"attributeSet":[{"attribute":{"value":"Sufficient evidence for dosage pathogenicity","dateValue":1341529200000,"type":"Haploinsufficiency"},"citation":[{"url":"http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=APC"}]},{"attribute":{"value":"No evidence available","dateValue":1341529200000,"type":"Triplosensitivity"},"citation":[{"url":"http://www.ncbi.nlm.nih.gov/projects/dbvar/ISCA/isca_gene.cgi?sym=APC"}]}],"sequenceLocation":[{"assembly":"GRCh38","chr":"5","accession":"NC_000005.10","start":112707504,"stop":112846238,"displayStart":112707504,"displayStop":112846238,"strand":"+","variantLength":138735,"assemblyAccessionVersion":"GCF_000001405.28","assemblyStatus":"current"},{"assembly":"GRCh37","chr":"5","accession":"NC_000005.9","start":112043201,"stop":112181935,"displayStart":112043201,"displayStop":112181935,"strand":"+","variantLength":138735,"assemblyAccessionVersion":"GCF_000001405.25","assemblyStatus":"previous"}],"comment":[{"value":"This gene is cited in the ACMG recommendations of 2013 (PubMed 23788249) for reporting incidental findings in exons.","dataSource":"NCBI curation","type":"PUBLIC"}],"type":"variant in gene","xref":[{"db":"Gene","id":"324","status":"CURRENT"},{"db":"OMIM","id":"611731","type":"MIM","status":"CURRENT"},{"db":"HGNC","id":"HGNC:583","status":"CURRENT"}]}],"type":"Deletion","id":15851,"xref":[{"db":"OMIM","id":"611731.0020","type":"Allelic variant","status":"CURRENT"},{"db":"dbSNP","id":"387906235","type":"rs","status":"CURRENT"}]}],"name":[{"elementValue":{"value":"NM_000038.5(APC):c.4391_4394delAGAG (p.Glu1464Valfs)","type":"Preferred"}}],"type":"Variant","id":812},"traitSet":{"trait":[{"name":[{"elementValue":{"value":"Periampullary adenoma","type":"Preferred"}},{"elementValue":{"value":"EFO:0000232","type":"EFO id"}},{"elementValue":{"value":"adenoma","type":"EFO name"}},{"elementValue":{"value":"http://www.ebi.ac.uk/efo/EFO_0000232","type":"EFO URL"}}],"attributeSet":[{"attribute":{"value":"Neoplasm","type":"keyword"}}],"type":"Disease","id":9669,"xref":[{"db":"MedGen","id":"CN068444","status":"CURRENT"}]}],"type":"Disease","id":210},"dateCreated":1344812400000,"dateLastUpdated":1455667200000,"id":58354},"clinVarAssertion":[{"clinVarSubmissionID":{"submitter":"OMIM","title":"APC, 4-BP DEL, CODON 1464_ADENOMA, PERIAMPULLARY, SOMATIC","localKey":"611731.0020_ADENOMA, PERIAMPULLARY, SOMATIC","submitterDate":1365030000000},"clinVarAccession":{"acc":"SCV000021001","version":1,"type":"SCV","orgID":3,"dateUpdated":1444518000000},"recordStatus":"current","clinicalSignificance":{"reviewStatus":"NO_ASSERTION_CRITERIA_PROVIDED","description":["Pathogenic"],"dateLastEvaluated":752112000000},"assertion":{"type":"variation to disease"},"externalID":{"db":"OMIM","id":"611731.0020","type":"Allelic variant","status":"CURRENT"},"observedIn":[{"sample":{"origin":"somatic","species":{"value":"human"},"affectedStatus":"not provided"},"method":[{"methodType":"LITERATURE_ONLY"}],"observedData":[{"attribute":{"value":"In tumor tissue of a periampullary adenoma from a patient with FAP (175100), Bapat et al. (1993) identified a somatic 4-bp deletion (AGAG) at codon 1464 of the APC gene. The patient had a germline APC mutation (611731.0023).","type":"Description"},"citation":[{"id":[{"value":"8281160","source":"PubMed"}]}],"xref":[{"db":"OMIM","id":"175100","type":"MIM","status":"CURRENT"}]}]}],"measureSet":{"measure":[{"name":[{"elementValue":{"value":"APC, 4-BP DEL, CODON 1464","type":"Preferred"}}],"attributeSet":[{"attribute":{"value":"4-BP DEL, CODON 1464","type":"NonHGVS"}}],"measureRelationship":[{"symbol":[{"elementValue":{"value":"APC","type":"Preferred"}}],"type":"variant in gene"}],"type":"Variation","xref":[{"db":"OMIM","id":"611731.0020","type":"Allelic variant","status":"CURRENT"}]}],"type":"Variant"},"traitSet":{"trait":[{"name":[{"elementValue":{"value":"ADENOMA, PERIAMPULLARY, SOMATIC","type":"Preferred"}}],"type":"Disease"}],"type":"Disease"},"id":21001}],"id":10311588})

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
        test_dict = {"literature": {"references": [{"lit_id": "http://europepmc.org/abstract/MED/8281160"}]}, "disease": {"id": ["http://www.ebi.ac.uk/efo/EFO_0000232"]}, "validated_against_schema_version": "1.2.3", "target": {"target_type": "http://identifiers.org/cttv.target/gene_variant", "id": ["http://identifiers.org/ensembl/ENSG00000134982"], "activity": "http://identifiers.org/cttv.activity/damaging_to_target"}, "sourceID": "eva_somatic", "type": "somatic_mutation", "access_level": "public", "unique_association_fields": {"gene": "ENSG00000134982", "alleleOrigin": "somatic", "phenotype": "http://www.ebi.ac.uk/efo/EFO_0000232", "clinvarAccession": "RCV000000851"}, "evidence": {"is_associated": True, "provenance_type": {"literature": {"references": [{"lit_id": "http://europepmc.org/abstract/MED/8281160"}]}, "expert": {"status": True, "statement": "Primary submitter of data"}, "database": {"id": "EVA", "dbxref": {"url": "http://identifiers.org/clinvar.record/RCV000000851", "id": "http://identifiers.org/clinvar", "version": "2015-04"}, "version": "1.0"}}, "evidence_codes": ["http://purl.obolibrary.org/obo/ECO_0000205"], "date_asserted": "2016-02-17T00:00:00", "urls": [{"url": "http://www.ncbi.nlm.nih.gov/clinvar/RCV000000851", "nice_name": "Further details in ClinVar database"}], "known_mutations": [{"preferred_name": "frameshift_variant", "functional_consequence": "http://purl.obolibrary.org/obo/SO_0001589"}], "resource_score": {"type": "probability", "value": 1}}}

        test_ev_string = evidence_strings.CTTVEvidenceString(test_dict)

        self.assertEqual(set(self.evidence_string.__dict__.values()), set(test_ev_string.__dict__.values()))


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
