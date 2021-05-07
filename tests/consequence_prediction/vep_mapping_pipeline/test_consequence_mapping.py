from consequence_prediction.vep_mapping_pipeline.consequence_mapping import extract_consequences


def test_extract_consequences():
    vep_results = [
        {
            'input': '10 27169969 . C A',
            'transcript_consequences': [
                {
                    'gene_symbol': 'MASTL',
                    'gene_id': 'ENSG00000120539',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['missense_variant']
                },
                # {  # missing gene symbol is okay
                #     'gene_id': 'ENSG00000120538',
                #     'biotype': 'protein_coding',
                #     'consequence_terms': ['missense_variant']
                # },
                {  # less severe consequences should be filtered out
                    'gene_symbol': 'ACBD5',
                    'gene_id': 'ENSG00000107897',
                    'biotype': 'protein_coding',
                    'consequence_terms': ['3_prime_UTR_variant']
                },
                {  # other biotypes should be filtered out
                    'gene_symbol': 'MIR125A',
                    'gene_id': 'ENSG00000208008',
                    'biotype': 'miRNA',
                    'consequence_terms': ['missense_variant']
                }
            ]
        }
    ]
    results = {}
    extract_consequences(
        vep_results=vep_results,
        acceptable_biotypes=['protein_coding'],
        only_closest=False,
        results_by_variant=results,
        report_distance=False,
    )
    assert results == {'10 27169969 . C A': [
        ('10 27169969 . C A', 'ENSG00000120539', 'MASTL', 'missense_variant', 0),
        # ('10 27169969 . C A', 'ENSG00000120538', '', 'missense_variant', 0),
    ]}


def test_extract_consequences_only_closest():
    vep_results = [
        {
            'input': '6 1611781 . C CACGGCG',
            'transcript_consequences': [
                {
                    'gene_symbol': 'FOXC1',
                    'gene_id': 'ENSG00000054598',
                    'biotype': 'protein_coding',
                    'distance': 42,
                    'consequence_terms': ['downstream_gene_variant'],
                },
                {  # more distant consequence of the same severity should be excluded
                    'gene_symbol': 'FOXCUT',
                    'gene_id': 'ENSG00000280916',
                    'biotype': 'protein_coding',
                    'distance': 4427,
                    'consequence_terms': ['downstream_gene_variant']
                }
            ],
        }
    ]
    results = {}
    extract_consequences(
        vep_results=vep_results,
        acceptable_biotypes=['protein_coding'],
        only_closest=True,
        results_by_variant=results,
        report_distance=True,
    )
    assert results == {'6 1611781 . C CACGGCG': [
        ('6 1611781 . C CACGGCG', 'ENSG00000054598', 'FOXC1', 'downstream_gene_variant', 42)
    ]}
