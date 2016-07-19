from datetime import datetime
import os
import unittest

from eva_cttv_pipeline import clinvar
from eva_cttv_pipeline import consequence_type as CT
from tests import test_clinvar_to_evidence_strings

import tests.test_config as test_config


class TestClinvarRecord(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_clinvar_record = get_test_record()
        cls.rcv_to_rs, cls.rcv_to_nsv = \
            clinvar.get_rcv_to_rsnsv_mapping(test_config.variant_summary_file)
        cls.consequence_type_dict = CT.process_consequence_type_file(test_config.snp_2_gene_file)

    def test_gene_id(self):
        self.assertEqual(self.test_clinvar_record.gene_id, "NM_152443")

    def test_ensembl_id(self):
        self.assertEqual(self.test_clinvar_record.ensembl_id, "ENSG00000072121")

    def test_date(self):
        self.assertEqual(self.test_clinvar_record.date,
                         datetime.fromtimestamp(1414627200000/1000).isoformat())

    def test_score(self):
        self.assertEqual(self.test_clinvar_record.score, 1)

    def test_acc(self):
        self.assertEqual(self.test_clinvar_record.accession, "RCV000002127")

    def test_traits(self):
        self.assertEqual(self.test_clinvar_record.traits, [['Leber congenital amaurosis 13',
                                                            'Leber Congenital Amaurosis']])

    def test_trait_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.trait_pubmed_refs, [[20301475, 20301590]])

    def test_observed_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.observed_pubmed_refs, [15258582, 15322982])

    def test_measure_set_pubmed_refs(self):
        self.assertEqual(self.test_clinvar_record.measure_set_pubmed_refs, [])

    def test_hgvs(self):
        self.assertEqual(self.test_clinvar_record.hgvs,
                         ['NM_152443.2:c.677A>G',
                          'NG_008321.1:g.32324A>G',
                          'NC_000014.9:g.67729209A>G',
                          'NC_000014.8:g.68195926A>G',
                          'NP_689656.2:p.Tyr226Cys'])

    def test_clinical_significance(self):
        self.assertEqual(self.test_clinvar_record.clinical_significance, "pathogenic")

    def test_get_rs(self):
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_rs(self.rcv_to_rs), "rs28940313")
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_rs({}), None)

    def test_get_nsv(self):
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_nsv(self.rcv_to_nsv), None)
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_nsv({"RCV000002127": "nsv123test"}),
                         "nsv123test")

    def test___get_main_consequence_types(self):
        test_consequence_type = CT.ConsequenceType(ensembl_gene_ids=["ENSG00000139988"],
                                                   so_names=["sequence_variant"])

        # print([so_name.__dict__ for so_name in self.consequence_type_dict["rs28940313"].so_terms])
        # print(self.rcv_to_rs)

        self.assertEqual(
            self.test_clinvar_record._ClinvarRecord__get_main_consequence_types(
                self.consequence_type_dict, self.rcv_to_rs),
            test_consequence_type)
        self.assertEqual(self.test_clinvar_record._ClinvarRecord__get_main_consequence_types(
            {}, {}),
            None)

    def test_variant_type(self):
        self.assertEqual(self.test_clinvar_record.variant_type, "single nucleotide variant")

    def test_allele_origins(self):
        self.assertEqual(self.test_clinvar_record.allele_origins, ['germline'])


class TestGetRcvToRSNSVMapping(unittest.TestCase):
    variant_summary_file_path = os.path.join(os.path.dirname(__file__), 'resources',
                                             'variant_summary_2015-05_test_extract.txt')
    rcv_to_rs, rcv_to_nsv = clinvar.get_rcv_to_rsnsv_mapping(variant_summary_file_path)

    def test_rcv_to_rs(self):
        self.assertEqual(self.rcv_to_rs["RCV000000012"], "rs397704705")
        self.assertEqual(self.rcv_to_rs["RCV000000381"], "rs137854556")
        self.assertEqual(self.rcv_to_rs["RCV000000204"], "rs121965059")

    def test_rcv_to_nsv(self):
        self.assertEqual(self.rcv_to_nsv["RCV000004182"], "nsv1067860")
        self.assertEqual(self.rcv_to_nsv["RCV000004183"], "nsv1067861")
        self.assertEqual(self.rcv_to_nsv["RCV000004554"], "nsv1067916")


def get_test_record():
    test_record = clinvar.ClinvarRecord(test_clinvar_to_evidence_strings.MAPPINGS,
        {
  "recordStatus": "current",
  "title": "NM_152443.2(RDH12):c.677A>G (p.Tyr226Cys) AND Leber congenital amaurosis 13",
  "clinVarAssertion": [
    {
      "measureSet": {
        "type": "Variant",
        "measure": [
          {
            "type": "Variation",
            "measureRelationship": [
              {
                "type": "variant in gene",
                "symbol": [
                  {
                    "elementValue": {
                      "type": "Preferred",
                      "value": "RDH12"
                    }
                  }
                ]
              }
            ],
            "name": [
              {
                "elementValue": {
                  "type": "Preferred",
                  "value": "RDH12, TYR226CYS"
                }
              }
            ],
            "xref": [
              {
                "db": "OMIM",
                "type": "Allelic variant",
                "status": "CURRENT",
                "id": "608830.0001"
              }
            ],
            "attributeSet": [
              {
                "attribute": {
                  "type": "NonHGVS",
                  "value": "TYR226CYS"
                }
              }
            ]
          }
        ]
      },
      "clinVarAccession": {
        "type": "SCV",
        "version": 1,
        "acc": "SCV000022285",
        "dateUpdated": 1414627200000,
        "orgID": 3
      },
      "assertion": {
        "type": "variation to disease"
      },
      "clinicalSignificance": {
        "dateLastEvaluated": 1354060800000,
        "description": [
          "Pathogenic"
        ]
      },
      "recordStatus": "current",
      "id": 22285,
      "observedIn": [
        {
          "sample": {
            "affectedStatus": "not provided",
            "species": {
              "value": "human"
            },
            "origin": "germline"
          },
          "observedData": [
            {
              "citation": [
                {
                  "id": {
                    "value": "15258582",
                    "source": "PubMed"
                  }
                }
              ],
              "attribute": {
                "type": "Description",
                "value": "In affected members of 3 consanguineous Austrian kindreds with Leber congenital amaurosis-13 (612712), Janecke et al. (2004) identified homozygosity for a 677A-G transition in exon 6 of the RDH12 gene, resulting in a tyr226-to-cys (Y226C) substitution. The same mutation was identified in 2 Austrian individuals with sporadic LCA13. Janecke et al. (2004) demonstrated that, when expressed in COS-7 cells, the cys226 variant had diminished activity in interconverting isomers of retinol and retinal."
              },
              "xref": [
                {
                  "db": "OMIM",
                  "type": "MIM",
                  "status": "CURRENT",
                  "id": "612712"
                }
              ]
            },
            {
              "citation": [
                {
                  "id": {
                    "value": "15322982",
                    "source": "PubMed"
                  }
                }
              ],
              "attribute": {
                "type": "Description",
                "value": "In affected members of a French family with LCA, Perrault et al. (2004) identified the Y226C mutation in compound heterozygous state with a 523T-C transition in exon 5 of the RDH12 gene, resulting in a ser175-to-pro substitution (S175P; 608830.0011)."
              }
            }
          ],
          "method": [
            {
              "methodType": "LITERATURE_ONLY"
            }
          ]
        }
      ],
      "traitSet": {
        "type": "Disease",
        "trait": [
          {
            "type": "Disease",
            "name": [
              {
                "elementValue": {
                  "type": "Preferred",
                  "value": "LEBER CONGENITAL AMAUROSIS 13"
                }
              }
            ]
          }
        ]
      },
      "clinVarSubmissionID": {
        "submitterDate": 1354060800000,
        "localKey": "608830.0001_LEBER CONGENITAL AMAUROSIS 13",
        "submitter": "OMIM",
        "title": "RDH12, TYR226CYS_LEBER CONGENITAL AMAUROSIS 13"
      },
      "externalID": {
        "db": "OMIM",
        "type": "Allelic variant",
        "status": "CURRENT",
        "id": "608830.0001"
      }
    }
  ],
  "referenceClinVarAssertion": {
    "dateCreated": 1344812400000,
    "measureSet": {
      "type": "Variant",
      "measure": [
        {
          "sequenceLocation": [
            {
              "start": 68195926,
              "alternateAllele": "G",
              "variantLength": 1,
              "chr": "14",
              "accession": "NC_000014.8",
              "stop": 68195926,
              "assembly": "GRCh37",
              "referenceAllele": "A"
            },
            {
              "start": 67729209,
              "alternateAllele": "G",
              "variantLength": 1,
              "chr": "14",
              "accession": "NC_000014.9",
              "stop": 67729209,
              "assembly": "GRCh38",
              "referenceAllele": "A"
            }
          ],
          "type": "single nucleotide variant",
          "name": [
            {
              "elementValue": {
                "type": "Preferred",
                "value": "NM_152443.2(RDH12):c.677A>G (p.Tyr226Cys)"
              }
            }
          ],
          "xref": [
            {
              "db": "OMIM",
              "type": "Allelic variant",
              "status": "CURRENT",
              "id": "608830.0001"
            },
            {
              "db": "dbSNP",
              "type": "rs",
              "status": "CURRENT",
              "id": "28940313"
            }
          ],
          "measureRelationship": [
            {
              "sequenceLocation": [
                {
                  "start": 68213236,
                  "chr": "14",
                  "accession": "NC_000014.8",
                  "stop": 68283305,
                  "assembly": "GRCh37",
                  "strand": "-"
                },
                {
                  "start": 67728891,
                  "chr": "14",
                  "accession": "NC_000014.9",
                  "stop": 67816589,
                  "assembly": "GRCh38",
                  "strand": "-"
                }
              ],
              "type": "variant in gene",
              "name": [
                {
                  "elementValue": {
                    "type": "Preferred",
                    "value": "zinc finger, FYVE domain containing 26"
                  }
                }
              ],
              "xref": [
                {
                  "db": "Gene",
                  "status": "CURRENT",
                  "id": "23503"
                },
                {
                  "db": "OMIM",
                  "type": "MIM",
                  "status": "CURRENT",
                  "id": "612012"
                }
              ],
              "symbol": [
                {
                  "elementValue": {
                    "type": "Preferred",
                    "value": "ZFYVE26"
                  }
                }
              ]
            },
            {
              "sequenceLocation": [
                {
                  "start": 68168602,
                  "chr": "14",
                  "accession": "NC_000014.8",
                  "stop": 68201167,
                  "assembly": "GRCh37",
                  "strand": "+"
                },
                {
                  "start": 67701885,
                  "chr": "14",
                  "accession": "NC_000014.9",
                  "stop": 67734450,
                  "assembly": "GRCh38",
                  "strand": "+"
                }
              ],
              "type": "variant in gene",
              "name": [
                {
                  "elementValue": {
                    "type": "Preferred",
                    "value": "retinol dehydrogenase 12 (all-trans/9-cis/11-cis)"
                  }
                }
              ],
              "xref": [
                {
                  "db": "Gene",
                  "status": "CURRENT",
                  "id": "145226"
                },
                {
                  "db": "OMIM",
                  "type": "MIM",
                  "status": "CURRENT",
                  "id": "608830"
                }
              ],
              "symbol": [
                {
                  "elementValue": {
                    "type": "Preferred",
                    "value": "RDH12"
                  }
                }
              ]
            }
          ],
          "cytogeneticLocation": [
            "14q24.1"
          ],
          "attributeSet": [
            {
              "attribute": {
                "type": "HGVS, coding, RefSeq",
                "value": "NM_152443.2:c.677A>G",
                "change": "c.677A>G"
              }
            },
            {
              "attribute": {
                "type": "HGVS, genomic, RefSeqGene",
                "value": "NG_008321.1:g.32324A>G",
                "change": "g.32324A>G"
              }
            },
            {
              "attribute": {
                "type": "HGVS, genomic, top level",
                "value": "NC_000014.9:g.67729209A>G",
                "integerValue": 38,
                "change": "g.67729209A>G"
              }
            },
            {
              "attribute": {
                "type": "HGVS, genomic, top level, previous",
                "value": "NC_000014.8:g.68195926A>G",
                "integerValue": 37,
                "change": "g.68195926A>G"
              }
            },
            {
              "attribute": {
                "type": "HGVS, protein, RefSeq",
                "value": "NP_689656.2:p.Tyr226Cys",
                "change": "p.Tyr226Cys"
              },
              "xref": [
                {
                  "db": "dbSNP",
                  "type": "rs",
                  "status": "CURRENT",
                  "id": "28940313"
                }
              ]
            },
            {
              "attribute": {
                "type": "MolecularConsequence",
                "value": "missense variant"
              },
              "xref": [
                {
                  "db": "Sequence Ontology",
                  "status": "CURRENT",
                  "id": "SO:0001583"
                },
                {
                  "db": "RefSeq",
                  "status": "CURRENT",
                  "id": "NM_152443.2:c.677A>G"
                }
              ]
            },
            {
              "attribute": {
                "type": "ProteinChange1LetterCode",
                "value": "Y226C"
              },
              "xref": [
                {
                  "db": "OMIM",
                  "type": "Allelic variant",
                  "status": "CURRENT",
                  "id": "608830.0001"
                }
              ]
            },
            {
              "attribute": {
                "type": "ProteinChange3LetterCode",
                "value": "TYR226CYS"
              },
              "xref": [
                {
                  "db": "OMIM",
                  "type": "Allelic variant",
                  "status": "CURRENT",
                  "id": "608830.0001"
                }
              ]
            }
          ],
          "id": 17085
        }
      ],
      "name": [
        {
          "elementValue": {
            "type": "preferred name",
            "value": "NM_152443.2(RDH12):c.677A>G (p.Tyr226Cys)"
          }
        }
      ],
      "id": 2046
    },
    "clinVarAccession": {
      "type": "RCV",
      "version": 1,
      "acc": "RCV000002127",
      "dateUpdated": 1414627200000
    },
    "assertion": {
      "type": "VARIATION_TO_DISEASE"
    },
    "clinicalSignificance": {
      "reviewStatus": "CLASSIFIED_BY_SINGLE_SUBMITTER",
      "dateLastEvaluated": 1354060800000,
      "description": "Pathogenic"
    },
    "recordStatus": "current",
    "id": 59630,
    "traitSet": {
      "type": "Disease",
      "id": 522,
      "trait": [
        {
          "type": "Disease",
          "name": [
            {
              "elementValue": {
                "type": "Preferred",
                "value": "Leber congenital amaurosis 13"
              },
              "xref": [
                {
                  "db": "Genetic Alliance",
                  "status": "CURRENT",
                  "id": "Leber+congenital+amaurosis+type+13/4135"
                }
              ]
            },
            {
              "elementValue": {
                "type": "Alternate",
                "value": "Leber Congenital Amaurosis"
              },
              "xref": [
                {
                  "db": "GeneReviews",
                  "status": "CURRENT",
                  "id": "NBK1298"
                }
              ]
            }
          ],
          "xref": [
            {
              "db": "MedGen",
              "status": "CURRENT",
              "id": "C2675186"
            },
            {
              "db": "OMIM",
              "type": "MIM",
              "status": "CURRENT",
              "id": "612712"
            }
          ],
          "attributeSet": [
            {
              "attribute": {
                "type": "public definition",
                "value": "Leber congenital amaurosis (LCA), a severe dystrophy of the retina, typically becomes evident in the first year of life. Visual function is usually poor and often accompanied by nystagmus, sluggish or near-absent pupillary responses, photophobia, high hyperopia, and keratoconus. Visual acuity is rarely better than 20/400. A characteristic finding is Franceschetti's oculo-digital sign, comprising eye poking, pressing, and rubbing. The appearance of the fundus is extremely variable. While the retina may initially appear normal, a pigmentary retinopathy reminiscent of retinitis pigmentosa is frequently observed later in childhood. The electroretinogram (ERG) is characteristically \"nondetectable\" or severely subnormal."
              },
              "xref": [
                {
                  "db": "GeneReviews",
                  "status": "CURRENT",
                  "id": "NBK1298"
                }
              ]
            }
          ],
          "id": 2542,
          "citation": [
            {
              "type": "review",
              "abbrev": "GeneReviews",
              "id": {
                "value": "20301475",
                "source": "PubMed"
              }
            },
            {
              "type": "review",
              "abbrev": "GeneReviews",
              "id": {
                "value": "20301590",
                "source": "PubMed"
              }
            }
          ],
          "symbol": [
            {
              "elementValue": {
                "type": "Preferred",
                "value": "LCA13"
              },
              "xref": [
                {
                  "db": "OMIM",
                  "type": "MIM",
                  "status": "CURRENT",
                  "id": "612712"
                }
              ]
            }
          ]
        }
      ]
    },
    "observedIn": [
      {
        "sample": {
          "affectedStatus": "not provided",
          "species": {
            "value": "human",
            "taxonomyId": 9606
          },
          "origin": "germline"
        },
        "observedData": [
          {
            "citation": [
              {
                "type": "general",
                "id": {
                  "value": "15258582",
                  "source": "PubMed"
                }
              }
            ],
            "attribute": {
              "type": "Description",
              "value": "In affected members of 3 consanguineous Austrian kindreds with Leber congenital amaurosis-13 (612712), Janecke et al. (2004) identified homozygosity for a 677A-G transition in exon 6 of the RDH12 gene, resulting in a tyr226-to-cys (Y226C) substitution. The same mutation was identified in 2 Austrian individuals with sporadic LCA13. Janecke et al. (2004) demonstrated that, when expressed in COS-7 cells, the cys226 variant had diminished activity in interconverting isomers of retinol and retinal."
            },
            "id": 3703574
          },
          {
            "citation": [
              {
                "type": "general",
                "id": {
                  "value": "15322982",
                  "source": "PubMed"
                }
              }
            ],
            "attribute": {
              "type": "Description",
              "value": "In affected members of a French family with LCA, Perrault et al. (2004) identified the Y226C mutation in compound heterozygous state with a 523T-C transition in exon 5 of the RDH12 gene, resulting in a ser175-to-pro substitution (S175P; 608830.0011)."
            },
            "id": 3703574
          }
        ],
        "method": [
          {
            "methodType": "LITERATURE_ONLY"
          }
        ]
      }
    ],
    "dateLastUpdated": 1414627200000
  },
  "id": 3908085
})
    # record_string = json.load(test_record)
    return test_record

